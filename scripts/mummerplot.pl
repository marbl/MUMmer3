#!__PERL_PATH

################################################################################
#   Programmer: Adam M Phillippy, The Institute for Genomic Research
#         File: mummerplot
#         Date: 01 / 08 / 03
#               01 / 06 / 05 rewritten (v3.0)
#  
#        Usage:
#    mummerplot  [options]  <match file>
# 
#                Try 'mummerplot -h' for more information.
# 
#      Purpose: To generate a gnuplot plot for the display of mummer, nucmer,
#               promer, and show-tiling alignments.
# 
################################################################################

use lib "__SCRIPT_DIR";
use Foundation;
use strict;


#================================================================= Globals ====#
#-- terminal types
my $X11    = "x11";
my $PS     = "postscript";
my $PNG    = "png";

#-- terminal sizes
my $SMALL  = "small";
my $MEDIUM = "medium";
my $LARGE  = "large";

my %TERMSIZE =
    (
     $X11 => { $SMALL => 500, $MEDIUM => 700,  $LARGE => 900  }, # screen pix
     $PS  => { $SMALL => 1,   $MEDIUM => 2,    $LARGE => 3    }, # pages
     $PNG => { $SMALL => 800, $MEDIUM => 1024, $LARGE => 1400 }  # image pix
     );

#-- terminal format
my $FFACE   = "Courier";
my $FSIZE   = "8";
my $FFORMAT = "%.0f";

#-- output suffixes
my $DATA    = "data";
my $GNUPLOT = "gnuplot";

my %SUFFIX =
    (
     $DATA    => ".plot",
     $GNUPLOT => ".gp",
     $PS      => ".ps",
     $PNG     => ".png"
     );


#================================================================= Options ====#
my $OPT_Mfile;                     # match file

my $OPT_coverage;                  # -c option
my $OPT_prefix   = "out";          # -p option
my $OPT_terminal = $X11;           # -t option
my $OPT_IdR;                       # -r option
my $OPT_IdQ;                       # -q option
my $OPT_Rfile;                     # -R option
my $OPT_Qfile;                     # -Q option
my $OPT_xrange;                    # -x option
my $OPT_yrange;                    # -y option

my $OPT_size     = $SMALL;         # -small, -medium, -large

my $OPT_Dfile;                     # .plot output
my $OPT_Gfile;                     # .gp output
my $OPT_Pfile;                     # .ps .png etc. output

my $OPT_gpstatus;                  # gnuplot status


#============================================================== Foundation ====#
my $VERSION = '3.0';

my $USAGE = qq~
  USAGE: mummerplot  [options]  <match file>
    ~;

my $HELP = qq~
  USAGE: mummerplot  [options]  <match file>

  DESCRIPTION:
    mummerplot generates plots of alignment data produced by mummer, nucmer,
    promer or show-tiling by using the GNU gnuplot utility. After generating
    the appropriate scripts and datafiles, mummerplot will attempt to run
    gnuplot to generate the plot. If this attempt fails, a warning will be
    output and the resulting .gp and .plot files will remain so that the user
    may run gnuplot independently. If the attempt succeeds, either an x11
    window will be spawned or an additional output file will be generated
    (.ps or .png depending on the selected terminal). Feel free to edit the
    resulting gnuplot files (.gp and .plot) and rerun gnuplot to change line
    thinkness, labels, colors, plot size etc.

  MANDATORY:
    match file      Set the alignment input to "match file"
                    Valid inputs are from mummer, nucmer, promer and
                    show-tiling (.out, .cluster, .delta and .tiling)

  OPTIONS:
    -c
    --[no]coverage  Generate a reference coverage plot (default for .tiling)
    --depend        Print the dependency information and exit
    -h
    --help          Display help information and exit
    -p|prefix       Set the prefix of the output files (default "$OPT_prefix")
    -r|IdR          Plot a particular reference sequence ID on the X-axis
    -q|IdQ          Plot a particular query sequence ID on the Y-axis
    -R|Rfile        Get the set of reference sequences to plot from Rfile
    -Q|Qfile        Get the set of query sequences to plot from Qfile
                    Rfile/Qfile Can either be the original DNA multi-FastA
                    files or lists of sequence IDs, lens and dirs [ /+/-]
    -s|size         Set the output terminal size to small, medium or large
                    (default "$OPT_size")
    -t|terminal     Set the output terminal to x11, postscript or png
                    (default "$OPT_terminal")
    -x|xrange       Set the xrange for the plot "[min:max]"
    -y|yrange       Set the yrange for the plot "[min:max]"
    -V
    --version       Display the version information and exit
    ~;

my @DEPEND = ("TIGR::Foundation", "gnuplot");

my $tigr = new TIGR::Foundation
    or die "ERROR: TIGR::Foundation could not be initialized\n";

$tigr -> setVersionInfo ($VERSION);
$tigr -> setUsageInfo ($USAGE);
$tigr -> setHelpInfo ($HELP);
$tigr -> addDependInfo (@DEPEND);


#=========================================================== Function Decs ====#
sub GetParseFunc( );

sub ParseIDs($$);

sub ParseDelta($);
sub ParseCluster($);
sub ParseMummer($);
sub ParseTiling($);

sub PlotData($$$);
sub WriteGP($$);
sub RunGP( );

sub ParseOptions( );


#=========================================================== Function Defs ====#
MAIN:
{
    my @aligns;                # (sR eR sQ eQ sim lenR lenQ idR idQ)
    my %refs;                  # (id => (off, len, [1/-1]))
    my %qrys;                  # (id => (off, len, [1/-1]))

    #-- Get the command line options (sets OPT_ global vars)
    ParseOptions( );

    #-- Get the alignment type
    my $parsefunc = GetParseFunc( );

    #-- Parse the reference and query IDs
    if    ( defined $OPT_IdR )   { $refs{$OPT_IdR} = [ 0, 0, 1 ]; }
    elsif ( defined $OPT_Rfile ) {
        ParseIDs ($OPT_Rfile, \%refs);
    }

    if    ( defined $OPT_IdQ )   { $qrys{$OPT_IdQ} = [ 0, 0, 1 ]; }
    elsif ( defined $OPT_Qfile ) {
        ParseIDs ($OPT_Qfile, \%qrys);
    }

    #-- Parse the alignment data
    $parsefunc->(\@aligns);
    
    #-- Plot the alignment data
    PlotData (\@aligns, \%refs, \%qrys);

    #-- Write the gnuplot script
    WriteGP (\%refs, \%qrys);

    #-- Run the gnuplot script
    unless ( $OPT_gpstatus == -1 ) {
        RunGP( );
    }

    exit (0);
}


#------------------------------------------------------------ GetParseFunc ----#
sub GetParseFunc ( )
{
    my $fref;

    open (MFILE, "<$OPT_Mfile")
        or die "ERROR: Could not open $OPT_Mfile, $!\n";

    $_ = <MFILE>;
    if ( !defined ) { die "ERROR: Could not read $OPT_Mfile, File is empty\n" }

  SWITCH: {
      #-- tiling
      if ( /^>\S+ \d+ bases/ ) {
          $fref = \&ParseTiling;
          last SWITCH;
      }

      #-- mummer
      if ( /^> \S+/ ) {
          $fref = \&ParseMummer;
          last SWITCH;
      }

      #-- nucmer/promer
      if ( /^\S+ \S+/ ) {
          $_ = <MFILE>;
          if ( (defined)  &&  (/^NUCMER$/  ||  /^PROMER$/) ) {
              $_ = <MFILE>;   # sequence header
              $_ = <MFILE>;   # alignment header
              if ( !defined ) {
                  $fref = \&ParseDelta;
                  last SWITCH;
              }
              elsif ( /^\d+ \d+ \d+ \d+ \d+ \d+ \d+$/ ) {
                  $fref = \&ParseDelta;
                  last SWITCH;
              }
              elsif ( /^[ \-][1-3] [ \-][1-3]$/ ) {
                  $fref = \&ParseCluster;
                  last SWITCH;
              }
          }
      }
      
      #-- default
      die "ERROR: Could not read $OPT_Mfile, Unrecognized file type\n";
  }

    close (MFILE)
        or print STDERR "WARNING: Trouble closing $OPT_Mfile, $!\n";

    return $fref;
}


#---------------------------------------------------------------- ParseIDs ----#
sub ParseIDs ($$)
{
    my $file = shift;
    my $href = shift;

    open (IDFILE, "<$file")
        or die "ERROR: Could not open $file, $!\n";

    my $dir;
    my $aref;
    my $isfasta;
    my $offset = 0;
    while ( <IDFILE> ) {
        #-- FastA header
        if ( /^>(\S+)/ ) {
            if ( exists $href->{$1} ) {
                print STDERR "WARNING: Duplicate sequence '$1' ignored\n";
                undef $aref;
                next;
            }

            if ( !$isfasta ) { $isfasta = 1; }
            if ( defined $aref ) { $offset += $aref->[1]; }

            $aref = [ $offset, 0, 1 ];
            $href->{$1} = $aref;
            next;
        }
        
        #-- FastA sequence
        if ( $isfasta  &&  /^\S+$/ ) {
            if ( defined $aref ) { $aref->[1] += (length) - 1; }
            next;
        }

        #-- ID len dir
        if ( !$isfasta  &&  /^(\S+)\s+(\d+)\s+([+-]?)$/ ) {
            if ( exists $href->{$1} ) {
                print STDERR "WARNING: Duplicate sequence '$1' ignored\n";
                undef $aref;
                next;
            }

            $dir = (defined $3 && $3 eq "-") ? -1 : 1;
            $aref = [ $offset, $2, $dir ];
            $offset += $2;
            $href->{$1} = $aref;
            next;
        }

        #-- default
        die "ERROR: Could not parse $file\n$_";
    }

    close (IDFILE)
        or print STDERR "WARNING: Trouble closing $file, $!\n";
}


#-------------------------------------------------------------- ParseDelta ----#
sub ParseDelta ($)
{
    my $aref = shift;

    open (MFILE, "<$OPT_Mfile")
        or die "ERROR: Could not open $OPT_Mfile, $!\n";

    print STDERR "Reading delta file $OPT_Mfile\n";

    my @align;
    my $ispromer;
    my ($sim, $tot);
    my ($lenR, $lenQ, $idR, $idQ);

    $_ = <MFILE>;
    $_ = <MFILE>;
    $ispromer = /^PROMER/;

    while ( <MFILE> ) {
        #-- delta int
        if ( /^([-]?\d+)$/ ) {
            if ( $1 < 0 ) {
                $tot ++;
            }
            elsif ( $1 == 0 ) {
                $align[4] = ($tot - $sim) / $tot * 100.0;
                push @$aref, [ @align ];
                $tot = 0;
            }
            next;
        }

        #-- alignment header
        if ( /^(\d+) (\d+) (\d+) (\d+) \d+ (\d+) \d+$/ ) {
            if ( $tot == 0 ) {
                @align = ($1, $2, $3, $4, 0, $lenR, $lenQ, $idR, $idQ);
                $tot = abs($1 - $2) + 1;
                $sim = $5;
                if ( $ispromer ) { $tot /= 3.0; }
                next;
            }
            #-- drop to default
        }

        #-- sequence header
        if ( /^>(\S+) (\S+) (\d+) (\d+)$/ ) {
            ($idR, $idQ, $lenR, $lenQ) = ($1, $2, $3, $4);
            $tot = 0;
            next;
        }

        #-- default
        die "ERROR: Could not parse $OPT_Mfile\n$_";
    }

    close (MFILE)
        or print STDERR "WARNING: Trouble closing $OPT_Mfile, $!\n";
}


#------------------------------------------------------------ ParseCluster ----#
sub ParseCluster ($)
{
    my $aref = shift;

    open (MFILE, "<$OPT_Mfile")
        or die "ERROR: Could not open $OPT_Mfile, $!\n";

    print STDERR "Reading cluster file $OPT_Mfile\n";

    my @align;
    my ($dR, $dQ, $len);
    my ($lenR, $lenQ, $idR, $idQ, $len);

    $_ = <MFILE>;
    $_ = <MFILE>;

    while ( <MFILE> ) {
        #-- match
        if ( /^\s+(\d+)\s+(\d+)\s+(\d+)\s+\S+\s+\S+$/ ) {
            @align = ($1, $1, $2, $2, 100, $lenR, $lenQ, $idR, $idQ);
            $len = $3 - 1;
            $align[1] += $dR == 1 ? $len : -$len;
            $align[3] += $dQ == 1 ? $len : -$len;
            push @$aref, [ @align ];
            next;
        }

        #-- cluster header
        if ( /^[ \-][1-3] [ \-][1-3]$/ ) {
            $dR = /^-/ ? -1 : 1;
            $dQ = /-[1-3]$/ ? -1 : 1;
            next;
        }

        #-- sequence header
        if ( /^>(\S+) (\S+) (\d+) (\d+)$/ ) {
            ($idR, $idQ, $lenR, $lenQ) = ($1, $2, $3, $4);
            next;
        }

        #-- default
        die "ERROR: Could not parse $OPT_Mfile\n$_";
    }

    close (MFILE)
        or print STDERR "WARNING: Trouble closing $OPT_Mfile, $!\n";
}


#------------------------------------------------------------- ParseMummer ----#
sub ParseMummer ($)
{
    my $aref = shift;

    open (MFILE, "<$OPT_Mfile")
        or die "ERROR: Could not open $OPT_Mfile, $!\n";

    print STDERR "Reading mummer file $OPT_Mfile (use mummer -c)\n";

    my @align;
    my ($dQ, $len);
    my ($lenQ, $idQ);

    while ( <MFILE> ) {
        #-- 3 column match
        if ( /^\s+(\d+)\s+(\d+)\s+(\d+)$/ ) {
            @align = ($1, $1, $2, $2, 100, 0, $lenQ, "REF", $idQ);
            $len = $3 - 1;
            $align[1] += $len;
            $align[3] += $dQ == 1 ? $len : -$len;
            push @$aref, [ @align ];
            next;
        }

        #-- 4 column match
        if ( /^\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)$/ ) {
            @align = ($2, $2, $3, $3, 100, 0, $lenQ, $1, $idQ);
            $len = $4 - 1;
            $align[1] += $len;
            $align[3] += $dQ == 1 ? $len : -$len;
            push @$aref, [ @align ];
            next;
        }

        #-- sequence header
        if ( /^> (\S+)/ ) {
            $idQ = $1;
            $dQ = /^> \S+ Reverse/ ? -1 : 1;
            $lenQ = /Len = (\d+)/ ? $1 : 0;
            next;
        }

        #-- default
        die "ERROR: Could not parse $OPT_Mfile\n$_";
    }

    close (MFILE)
        or print STDERR "WARNING: Trouble closing $OPT_Mfile, $!\n";
}


#------------------------------------------------------------- ParseTiling ----#
sub ParseTiling ($)
{
    my $aref = shift;

    if ( ! defined $OPT_coverage ) { $OPT_coverage = 1; }

    open (MFILE, "<$OPT_Mfile")
        or die "ERROR: Could not open $OPT_Mfile, $!\n";

    print STDERR "Reading tiling file $OPT_Mfile\n";

    my @align;
    my ($dR, $dQ, $len);
    my ($lenR, $lenQ, $idR, $idQ, $len);

    while ( <MFILE> ) {
        #-- tile
        if ( /^(\S+)\s+\S+\s+\S+\s+(\d+)\s+\S+\s+(\S+)\s+([+-])\s+(\S+)$/ ) {
            @align = ($1, $1, 1, 1, $3, $lenR, $2, $idR, $5);
            $len = $2 - 1;
            $align[1] += $len;
            $align[($4 eq "-" ? 2 : 3)] += $len;
            push @$aref, [ @align ];
            next;
        }

        #-- sequence header
        if ( /^>(\S+) (\d+) bases$/ ) {
            ($idR, $lenR) = ($1, $2);
            next;
        }

        #-- default
        die "ERROR: Could not parse $OPT_Mfile\n$_";
    }

    close (MFILE)
        or print STDERR "WARNING: Trouble closing $OPT_Mfile, $!\n";
}


#---------------------------------------------------------------- PlotData ----#
sub PlotData ($$$)
{
    my $aref = shift;
    my $rref = shift;
    my $qref = shift;

    open (DFILE, ">$OPT_Dfile")
        or die "ERROR: Could not open $OPT_Dfile, $!\n";

    my $align;
    my $isplotted;
    my $ismultiref;
    my $ismultiqry;
    my ($plenR, $plenQ, $pidR, $pidQ);
    my $index;

    foreach $index (0,1) {
        #-- index 0 == forward, index 1 == reverse
        print DFILE "#-- index $index\n";
        my $usedindex;

        foreach $align (@$aref) {

            my ($sR, $eR, $sQ, $eQ, $sim, $lenR, $lenQ, $idR, $idQ) = @$align;

            if ( ! defined $pidR ) {
                ($plenR, $plenQ, $pidR, $pidQ) = ($lenR, $lenQ, $idR, $idQ);
            }

            #-- set the sequence offset, length, direction, etc...
            my ($refoff, $reflen, $refdir);
            my ($qryoff, $qrylen, $qrydir);

            if ( defined (%$rref) ) {
                #-- skip reference sequence or set atts from hash
                if ( !exists ($rref->{$idR}) ) { next; }
                else { ($refoff, $reflen, $refdir) = @{$rref->{$idR}}; }
            }
            else {
                #-- no reference hash, so default atts
                ($refoff, $reflen, $refdir) = (0, 0, 1);
            }

            if ( defined (%$qref) ) {
                #-- skip query sequence or set atts from hash
                if ( !exists ($qref->{$idQ}) ) { next; }
                else { ($qryoff, $qrylen, $qrydir) = @{$qref->{$idQ}}; }
            }
            else {
                #-- no query hash, so default atts
                ($qryoff, $qrylen, $qrydir) = (0, 0, 1);
            }

            #-- get the coordinates right and plot the data
            if ( $refdir == -1 ) {
                $sR = $reflen - $sR + 1;
                $eR = $reflen - $eR + 1;
            }
            if ( $qrydir == -1 ) {
                $sQ = $qrylen - $sQ + 1;
                $eQ = $qrylen - $eQ + 1;
            }

            #-- skip for now if wrong index, 0 == forward, 1 == reverse
            unless ( $index == (($sR < $eR) != ($sQ < $eQ)) ) {
                next;
            }

            if ( $OPT_coverage ) {
                print DFILE ($sR + $refoff), " 10\n";
                print DFILE ($eR + $refoff), " 10\n\n";
                print DFILE ($sR + $refoff), " ", ($sim), "\n";
                print DFILE ($eR + $refoff), " ", ($sim), "\n\n";
            }
            else {
                print DFILE ($sR + $refoff), " ", ($sQ + $qryoff), "\n";
                print DFILE ($eR + $refoff), " ", ($eQ + $qryoff), "\n\n";
            }            

            #-- set some flags
            if ( !$ismultiref && $idR ne $pidR ) { $ismultiref = 1; }
            if ( !$ismultiqry && $idQ ne $pidQ ) { $ismultiqry = 1; }
            if ( !$isplotted ) { $isplotted = 1; }
            if ( !$usedindex ) { $usedindex = 1; }
        }
        #-- new index
        if ( !$usedindex ) {
            print DFILE "1 1\n1 1\n\n";
        }
        print DFILE "\n";
    }

    if ( !defined (%$rref) ) {
        if ( $ismultiref ) {
            print STDERR
                "WARNING: Multiple ref sequences overlaid, try -R or -r\n";
        }
        elsif ( defined $pidR ) {
            $rref->{$pidR} = [ 0, $plenR, 1 ];
        }
    }

    if ( !defined (%$qref) ) {
        if ( $ismultiqry && !$OPT_coverage ) {
            print STDERR
                "WARNING: Multiple qry sequences overlaid, try -Q, -q or -c\n";
        }
        elsif ( defined $pidQ ) {
            $qref->{$pidQ} = [ 0, $plenQ, 1 ];
        }
    }

    if ( !$isplotted ) {
        print STDERR "WARNING: No alignment data was plotted\n";
    }

    close (DFILE)
        or print STDERR "WARNING: Trouble closing $OPT_Dfile, $!\n";
}


#----------------------------------------------------------------- WriteGP ----#
sub WriteGP ($$)
{
    my $rref = shift;
    my $qref = shift;

    my $SIZE = $TERMSIZE{$OPT_terminal}{$OPT_size};

    my %TERMINAL = $OPT_gpstatus ?
        (
         $X11 => "set terminal $X11\n",
         $PS  => "set terminal $PS color solid \"$FFACE\" $FSIZE\n",
         $PNG => "set terminal $PNG small\n"
         )
        :
        (
         $X11 => "set terminal $X11 font \"$FFACE,$FSIZE\"\n",
         $PS  => "set terminal $PS color solid \"$FFACE\" $FSIZE\n",
         $PNG => "set terminal $PNG tiny size $SIZE,$SIZE\n"
         );
    
    my %PLOT = $OPT_coverage ?
        (
         $X11 => [ "w l lt 1 lw 3", "w l lt 3 lw 3" ],
         $PS  => [ "w l lt 1 lw 4", "w l lt 3 lw 4" ],
         $PNG => [ "w l lt 1 lw 3", "w l lt 3 lw 3" ]
         )
        :
        (
         $X11 => [ "w lp lt 1 lw 2 pt 6 ps 1.0", "w lp lt 3 lw 2 pt 6 ps 1.0" ],
         $PS  => [ "w lp lt 1 lw 2 pt 6 ps 0.5", "w lp lt 3 lw 2 pt 6 ps 0.5" ],
         $PNG => [ "w lp lt 1 lw 3 pt 6 ps 1.0", "w lp lt 3 lw 3 pt 6 ps 1.0" ]
         );

    open (GFILE, ">$OPT_Gfile")
        or die "ERROR: Could not open $OPT_Gfile, $!\n";

    my @refk = keys (%$rref);
    my @qryk = keys (%$qref);
    my ($xrange, $yrange);
    my ($xlabel, $ylabel);
    my ($tic, $dir);
    my $border = 0;

    #-- terminal header and output
    print GFILE $TERMINAL{$OPT_terminal};

    if ( defined $OPT_Pfile ) {
        print GFILE "set output \"$OPT_Pfile\"\n";
    }

    #-- set tics, determine labels, ranges (ref)
    if ( scalar (@refk) == 1 ) {
        $xlabel = $refk[0];
        $xrange = $rref->{$xlabel}[1];
    }
    else {
        $xrange = 0;
        print GFILE "set xtics rotate \( \\\n";
        foreach $xlabel ( sort { $rref->{$a}[0] <=> $rref->{$b}[0] } @refk ) {
            $xrange += $rref->{$xlabel}[1];
            $tic = $rref->{$xlabel}[0];
            $dir = ($rref->{$xlabel}[2] == 1) ? "" : "*";
            print GFILE "\"$dir$xlabel\" $tic, \\\n";
        }
        print GFILE "\"\" $xrange \\\n\)\n";
        $xlabel = "REF";
    }
    if ( $xrange == 0 ) { $xrange = "*"; }

    #-- set tics, determine labels, ranges (qry)
    if ( $OPT_coverage ) {
        $ylabel = "%SIM";
        $yrange = 110;
    }
    elsif ( scalar (@qryk) == 1 ) {
        $ylabel = $qryk[0];
        $yrange = $qref->{$ylabel}[1];
    }
    else {
        $yrange = 0;
        print GFILE "set ytics \( \\\n";
        foreach $ylabel ( sort { $qref->{$a}[0] <=> $qref->{$b}[0] } @qryk ) {
            $yrange += $qref->{$ylabel}[1];
            $tic = $qref->{$ylabel}[0];
            $dir = ($qref->{$ylabel}[2] == 1) ? "" : "*";
            print GFILE "\"$dir$ylabel\" $tic, \\\n";
        }
        print GFILE "\"\" $yrange \\\n\)\n";
        $ylabel = "QRY";
    }
    if ( $yrange == 0 ) { $yrange = "*"; }

    #-- determine borders
    if ( $xrange ne "*" && scalar (@refk) == 1 ) { $border |= 10; }
    if ( $yrange ne "*" && scalar (@qryk) == 1 ) { $border |= 5; }
    if ( $OPT_coverage ) { $border |= 5; }

    #-- grid, labels, border
    print GFILE
        "set grid\n",
        "set nokey\n",
        "set border $border\n",
        "set ticscale 0 0.5\n",
        "set xlabel \"$xlabel\"\n",
        "set ylabel \"$ylabel\"\n";

    #-- ranges
    if ( defined $OPT_xrange ) {
        print GFILE "set xrange $OPT_xrange\n";
    }
    else {
        print GFILE "set xrange [1:$xrange]\n";
    }

    if ( defined $OPT_yrange ) {
        print GFILE "set yrange $OPT_yrange\n";
    }
    else {
        print GFILE "set yrange [1:$yrange]\n";
    }

    #-- coords format
    print GFILE "set format \"$FFORMAT\"\n";
    if ( $OPT_gpstatus == 0 ) {
        print GFILE "set mouse format \"$FFORMAT\"\n";
    }

    #-- terminal specific sizes
    if ( $OPT_terminal eq $X11 ) {
        print GFILE ($OPT_coverage ? "set size 1,1\n" : "set size 1,1\n");
    }
    elsif ( $OPT_terminal eq $PS ) {
        if ( $OPT_coverage ) {
            print GFILE "set size ", 1.0000 * $SIZE, ",", 0.5 * $SIZE, "\n";
        }
        else {
            print GFILE "set size ", 0.7727 * $SIZE, ",", 1.0 * $SIZE, "\n";
        }
    }
    elsif ( $OPT_terminal eq $PNG ) {
        print GFILE ($OPT_coverage ? "set size 1,.375\n" : "set size 1,1\n");
    }

    #-- plot command
    print GFILE "plot \\\n",
    " \"$OPT_Dfile\" index 0 title \"fwd\" ${PLOT{$OPT_terminal}[0]}, \\\n",
    " \"$OPT_Dfile\" index 1 title \"rev\" ${PLOT{$OPT_terminal}[1]}\n";

    #-- interactive mode
    if ( $OPT_terminal eq $X11 ) {
        print GFILE "print \"-- INTERACTIVE MODE --\"\n";
        print GFILE "print \"consult gnuplot docs for command list\"\n";
        print GFILE "print \"left mouse for coords, right mouse to zoom\"\n";
        print GFILE "print \"press enter to exit\"\n";
        print GFILE "pause -1\n";
    }

    close (GFILE)
        or print STDERR "WARNING: Trouble closing $OPT_Gfile, $!\n";
}


#------------------------------------------------------------------- RunGP ----#
sub RunGP ( )
{
    if ( defined $OPT_Pfile ) {
        print STDERR "Rendering plot $OPT_Pfile\n";
    }
    else {
        print STDERR "Rendering plot to screen\n";
    }

    my $cmd = "gnuplot";
    if ( $OPT_terminal eq $X11 ) {
        my $size = $TERMSIZE{$OPT_terminal}{$OPT_size};
        $cmd .= " -geometry ${size}x";
        if ( $OPT_coverage ) { $size = sprintf ("%.0f", $size * .375); }
        $cmd .= "${size}+0+0 -title mummerplot";
    }
    $cmd .= " $OPT_Gfile";

    system ($cmd) and die "ERROR: Unable to run '$cmd', Unknown error\n";
}


#------------------------------------------------------------ ParseOptions ----#
sub ParseOptions ( )
{
    my ($opt_small, $opt_medium, $opt_large);
    my ($opt_ps, $opt_x11, $opt_png);
    my $opt_color;

    #-- Get options
    my $err = $tigr -> TIGR_GetOptions
        (
         "c|coverage!"  => \$OPT_coverage,
         "p|prefix=s"   => \$OPT_prefix,
         "r|IdR=s"      => \$OPT_IdR,
         "q|IdQ=s"      => \$OPT_IdQ,
         "R|Rfile=s"    => \$OPT_Rfile,
         "Q|Qfile=s"    => \$OPT_Qfile,
         "s|size=s"     => \$OPT_size,
         "t|terminal=s" => \$OPT_terminal,
         "x|xrange=s"   => \$OPT_xrange,
         "y|yrange=s"   => \$OPT_yrange,

         #-- hidden options
         "x11"          => \$opt_x11,
         "postscript"   => \$opt_ps,
         "png"          => \$opt_png,
         "small"        => \$opt_small,
         "medium"       => \$opt_medium,
         "large"        => \$opt_large,
         "color!"       => \$opt_color     # unsupported
         );

    if ( !$err  ||  scalar (@ARGV) != 1 ) {
        $tigr -> printUsageInfo( );
        die "Try '$0 -h' for more information.\n";
    }

    #-- Hidden options
    if ( defined $opt_color ) {
        print STDERR "WARNING: --[no]color option is no longer supported\n";
    }

    if    ( $opt_x11 ) { $OPT_terminal = $X11; }
    elsif ( $opt_ps  ) { $OPT_terminal = $PS;  }
    elsif ( $opt_png ) { $OPT_terminal = $PNG; }

    if    ( $opt_large  ) { $OPT_size = $LARGE;  }
    elsif ( $opt_medium ) { $OPT_size = $MEDIUM; }
    elsif ( $opt_small  ) { $OPT_size = $SMALL; }

    #-- Check options
    if ( !exists $TERMSIZE{$OPT_terminal} ) {
        die "ERROR: Invalid terminal type, $OPT_terminal\n";
    }

    if ( !exists $TERMSIZE{$OPT_terminal}{$OPT_size} ) {
        die "ERROR: Invalid terminal size, $OPT_size\n";
    }
    
    if ( $OPT_xrange ) {
        $OPT_xrange =~ tr/,/:/;
        $OPT_xrange =~ /^\[\d+:\d+\]$/
            or die "ERROR: Invalid xrange format, $OPT_xrange\n";
    }

    if ( $OPT_yrange ) {
        $OPT_yrange =~ tr/,/:/;
        $OPT_yrange =~ /^\[\d+:\d+\]$/
            or die "ERROR: Invalid yrange format, $OPT_yrange\n";
    }

    #-- Set file names
    $OPT_Mfile = $ARGV[0];
    $tigr->isReadableFile ($OPT_Mfile)
        or die "ERROR: Could not read $OPT_Mfile, $!\n";

    $OPT_Dfile = $OPT_prefix . $SUFFIX{$DATA};
    $tigr->isWritableFile ($OPT_Dfile) or $tigr->isCreatableFile ($OPT_Dfile)
        or die "ERROR: Could not write $OPT_Dfile, $!\n";

    $OPT_Gfile = $OPT_prefix . $SUFFIX{$GNUPLOT};
    $tigr->isWritableFile ($OPT_Gfile) or $tigr->isCreatableFile ($OPT_Gfile)
        or die "ERROR: Could not write $OPT_Gfile, $!\n";

    if ( exists $SUFFIX{$OPT_terminal} ) {
        $OPT_Pfile = $OPT_prefix . $SUFFIX{$OPT_terminal};
        $tigr->isWritableFile($OPT_Pfile) or $tigr->isCreatableFile($OPT_Pfile)
            or die "ERROR: Could not write $OPT_Pfile, $!\n";
    }

    !$OPT_Rfile or $tigr->isReadableFile ($OPT_Rfile)
        or die "ERROR: Could not read $OPT_Rfile, $!\n";

    !$OPT_Qfile or $tigr->isReadableFile ($OPT_Qfile)
        or die "ERROR: Could not read $OPT_Qfile, $!\n";

    #-- Check that status of gnuplot
    $OPT_gpstatus = system ("gnuplot --version");

    if ( $OPT_gpstatus == -1 ) {
        print STDERR
            "WARNING: Could not find gnuplot, plot will not be rendered\n";
    }
    elsif ( $OPT_gpstatus ) {
        print STDERR
            "WARNING: Using outdated gnuplot, use v4.0 for best results\n";
    }
}
