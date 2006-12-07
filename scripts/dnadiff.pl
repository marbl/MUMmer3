#!__PERL_PATH -w

#-------------------------------------------------------------------------------
#   Programmer: Adam M Phillippy, University of Maryland
#         File: dnadiff
#         Date: 11 / 29 / 06
#
#        Usage:
#    dnadiff  [options]  <Reference>  <Query>
#
#                Try 'dnadiff -h' for more information.
#
#      Purpose: Run comparative analysis of Reference vs. Query.
#
#-------------------------------------------------------------------------------

use lib "__SCRIPT_DIR";
use Foundation;
use File::Spec::Functions;
use strict;

my $BIN_DIR = "__BIN_DIR";
my $SCRIPT_DIR = "__SCRIPT_DIR";

my $VERSION_INFO = q~
DNAdiff version 1.00
    ~;

my $HELP_INFO = q~
  USAGE: dnadiff  [options]  <Reference>  <Query>
    or   dnadiff  [options]  -d <Delta File>

  DESCRIPTION:
    Run comparative analysis of two genomes using nucmer.

  MANDATORY:
    Reference       Set the input reference multi-FASTA filename
    Query           Set the input query multi-FASTA filename

  OPTIONS:
    -d|delta        Provide precomputed delta file for analysis
    -h
    --help          Display help information and exit
    -p|prefix       Set the prefix of the output files (default "out")
    -V
    --version       Display the version information and exit
    ~;


my $USAGE_INFO = q~
  USAGE: dnadiff  [options]  <Reference>  <Query>
    or   dnadiff  [options]  -d <Delta File>
    ~;


my @DEPEND_INFO =
    (
     "$BIN_DIR/delta-filter",
     "$BIN_DIR/show-diff",
     "$BIN_DIR/show-snps",
     "$BIN_DIR/show-coords",
     "$BIN_DIR/nucmer",
     "$SCRIPT_DIR/Foundation.pm"
     );

my $DELTA_FILTER = "$BIN_DIR/delta-filter";
my $SHOW_DIFF = "$BIN_DIR/show-diff";
my $SHOW_SNPS = "$BIN_DIR/show-snps";
my $SHOW_COORDS = "$BIN_DIR/show-coords";
my $NUCMER = "$BIN_DIR/nucmer";

my $SNPBuff         = 20;            # required buffer around "good" snps
my $OPT_Prefix      = "out";         # prefix for all output files
my $OPT_RefFile;                     # reference file
my $OPT_QryFile;                     # query file
my $OPT_DeltaFile;                   # unfiltered alignment file
my $OPT_DeltaFileO  = ".odelta";     # 1-to-1 delta alignment
my $OPT_DeltaFileM  = ".mdelta";     # M-to-M delta alignment
my $OPT_CoordsFileO = ".ocoords";    # 1-to-1 alignment coords
my $OPT_CoordsFileM = ".mcoords";    # M-to-M alignment coords
my $OPT_SnpsFile    = ".snps";       # snps output file
my $OPT_DiffRFile   = ".rdiff";      # diffile for R
my $OPT_DiffQFile   = ".qdiff";      # diffile for Q
my $OPT_ReportFile  = ".report";     # report file

my $TIGR;  # TIGR Foundation object


sub RunAlignment();
sub RunFilter();
sub RunCoords();
sub RunSNPs();
sub RunDiff();
sub MakeReport();

sub FastaSizes($$);

sub FileOpen($$);
sub FileClose($$);

sub GetOpt();


#--------------------------------------------------------------------- main ----
 main:
{
    GetOpt();

    RunAlignment() unless defined($OPT_DeltaFile);
    RunFilter();
    RunCoords();
    RunSNPs();
    RunDiff();
    MakeReport();

    exit(0);
}


#------------------------------------------------------------- RunAlignment ----
# Run nucmer
sub RunAlignment()
{
    print STDERR "Building alignments\n";
    my $cmd = "$NUCMER --maxmatch -p $OPT_Prefix $OPT_RefFile $OPT_QryFile";
    my $err = "ERROR: Failed to run nucmer, aborting.\n";

    system($cmd) == 0 or die $err;
    $OPT_DeltaFile = $OPT_Prefix . ".delta";
}


#---------------------------------------------------------------- RunFilter ----
# Run delta-filter
sub RunFilter()
{
    print STDERR "Filtering alignments\n";
    my $cmd1 = "$DELTA_FILTER -1 $OPT_DeltaFile > $OPT_DeltaFileO";
    my $cmd2 = "$DELTA_FILTER -m $OPT_DeltaFile > $OPT_DeltaFileM";
    my $err = "ERROR: Failed to run delta-filter, aborting.\n";

    system($cmd1) == 0 or die $err;
    system($cmd2) == 0 or die $err;
}


#------------------------------------------------------------------ RunSNPs ----
# Run show-snps
sub RunSNPs()
{
    print STDERR "Analyzing SNPs\n";
    my $cmd = "$SHOW_SNPS -rlTHC $OPT_DeltaFileO > $OPT_SnpsFile";
    my $err = "ERROR: Failed to run show-snps, aborting.\n";

    system($cmd) == 0 or die $err;
}


#---------------------------------------------------------------- RunCoords ----
# Run show-coords
sub RunCoords()
{
    print STDERR "Extracting alignment coordinates\n";
    my $cmd1 = "$SHOW_COORDS -rclTH $OPT_DeltaFileO > $OPT_CoordsFileO";
    my $cmd2 = "$SHOW_COORDS -rclTH $OPT_DeltaFileM > $OPT_CoordsFileM";
    my $err = "ERROR: Failed to run show-coords, aborting.\n";

    system($cmd1) == 0 or die $err;
    system($cmd2) == 0 or die $err;
}


#------------------------------------------------------------------ RunDiff ----
# Run show-diff
sub RunDiff()
{
    print STDERR "Extracting alignment breakpoints\n";
    my $cmd1 = "$SHOW_DIFF -rH $OPT_DeltaFileM > $OPT_DiffRFile";
    my $cmd2 = "$SHOW_DIFF -qH $OPT_DeltaFileM > $OPT_DiffQFile";
    my $err = "ERROR: Failed to run show-diff, aborting.\n";

    system($cmd1) == 0 or die $err;
    system($cmd2) == 0 or die $err;
}


#--------------------------------------------------------------- MakeReport ----
# Output alignment report
sub MakeReport()
{
    print STDERR "Generating report file\n";

    my ($fhi, $fho);
    my ($lo, $hi);
    my (%refs, %qrys) = ((),());
    my ($rqnAlignsO, $rqnAlignsM) = (0,0);
    my ($rSumLenO, $qSumLenO) = (0,0);
    my ($rSumLenM, $qSumLenM) = (0,0);
    my ($rqSumIdyO, $rqSumIdyM) = (0,0);
    my ($qnIns, $rnIns) = (0,0);
    my ($qSumIns, $rSumIns) = (0,0);
    my ($qnTIns, $rnTIns) = (0,0);
    my ($qSumTIns, $rSumTIns) = (0,0);
    my ($qnInv, $rnInv) = (0,0);
    my ($qnRel, $rnRel) = (0,0);
    my ($qnTrn, $rnTrn) = (0,0);
    my ($rnSeqs, $qnSeqs) = (0,0);
    my ($rnASeqs, $qnASeqs) = (0,0);
    my ($rnBases, $qnBases) = (0,0);
    my ($rnABases, $qnABases) = (0,0);
    my ($rnBrk, $qnBrk) = (0,0);
    my ($rqnSNPs, $rqnIndels) = (0,0);
    my ($rqnGSNPs, $rqnGIndels) = (0,0);
    my %rqSNPs =
        ( "."=>{"A"=>0,"C"=>0,"G"=>0,"T"=>0},
          "A"=>{"."=>0,"C"=>0,"G"=>0,"T"=>0},
          "C"=>{"."=>0,"A"=>0,"G"=>0,"T"=>0},
          "G"=>{"."=>0,"A"=>0,"C"=>0,"T"=>0},
          "T"=>{"."=>0,"A"=>0,"C"=>0,"G"=>0} );

    my %rqGSNPs =
        ( "."=>{"A"=>0,"C"=>0,"G"=>0,"T"=>0},
          "A"=>{"."=>0,"C"=>0,"G"=>0,"T"=>0},
          "C"=>{"."=>0,"A"=>0,"G"=>0,"T"=>0},
          "G"=>{"."=>0,"A"=>0,"C"=>0,"T"=>0},
          "T"=>{"."=>0,"A"=>0,"C"=>0,"G"=>0} );
    my $header;

    #-- Get delta header
    $fhi = FileOpen("<", $OPT_DeltaFile);
    $header .= <$fhi>;
    $header .= <$fhi>;
    $header .= "\n";
    FileClose($fhi, $OPT_DeltaFile);

    #-- Collect all reference and query IDs and lengths
    FastaSizes($OPT_RefFile, \%refs);
    FastaSizes($OPT_QryFile, \%qrys);

    #-- Count ref and qry seqs and lengths
    foreach ( values(%refs) ) {
        $rnSeqs++;
        $rnBases += $_;
    }
    foreach ( values(%qrys) ) {
        $qnSeqs++;
        $qnBases += $_;
    }

    #-- Count aligned seqs, aligned bases, and breakpoints for each R and Q
    $fhi = FileOpen("<", $OPT_CoordsFileM);
    while (<$fhi>) {
        chomp;
        my @A = split "\t";
        scalar(@A) == 13
            or die "ERROR: Unrecognized format $OPT_CoordsFileM, aborting.\n";

        $rqnAlignsM++;
        $rSumLenM += $A[4];
        $qSumLenM += $A[5];
        $rqSumIdyM += $A[6];

        if ( $refs{$A[11]} > 0 ) {
            $rnASeqs++;
            $rnABases += $refs{$A[11]};
            $refs{$A[11]} *= -1;
        }
        if ( $qrys{$A[12]} > 0 ) {
            $qnASeqs++;
            $qnABases += $qrys{$A[12]};
            $qrys{$A[12]} *= -1;
        }

        if ( $A[1] < $A[2] ) { $lo = $A[1]; $hi = $A[2]; }
        else                 { $lo = $A[2]; $hi = $A[1]; }
        $rnBrk++ if ( $lo != 1 );
        $rnBrk++ if ( $hi != $A[7] );

        if ( $A[3] < $A[4] ) { $lo = $A[3]; $hi = $A[4]; }
        else                 { $lo = $A[4]; $hi = $A[3]; }
        $qnBrk++ if ( $lo != 1 );
        $qnBrk++ if ( $hi != $A[8] );
    }
    FileClose($fhi, $OPT_CoordsFileM);

    #-- Calculate average %idy, length, etc.
    $fhi = FileOpen("<", $OPT_CoordsFileO);
    while (<$fhi>) {
        chomp;
        my @A = split "\t";
        scalar(@A) == 13
            or die "ERROR: Unrecognized format $OPT_CoordsFileO, aborting.\n";

        $rqnAlignsO++;
        $rSumLenO += $A[4];
        $qSumLenO += $A[5];
        $rqSumIdyO += $A[6];
    }
    FileClose($fhi, $OPT_CoordsFileO);

    #-- If you are reading this, you need to get out more...

    #-- Count reference diff features and indels
    $fhi = FileOpen("<", $OPT_DiffRFile);
    while (<$fhi>) {
        chomp;
        my @A = split "\t";
        defined($A[4])
            or die "ERROR: Unrecognized format $OPT_DiffRFile, aborting.\n";
        my $gap = $A[4];
        my $ins = $gap;

        #-- Check for tandem signature
        if ( $A[1] eq "GAP" ) {
            scalar(@A) == 7
                or die "ERROR: Unrecognized format $OPT_DiffRFile, aborting.\n";
            $ins = $A[6] if ( $A[6] > $gap );
            if ( $A[4] <= 0 && $A[5] <= 0 && $A[6] > 0 ) {
                $rnTIns++;
                $rSumTIns += $A[6];
             }
        }

        $rnABases -= $gap if ( $gap > 0 );
        $rnInv++ if ( $A[1] eq "INV" );
        $rnRel++ if ( $A[1] eq "JMP" );
        $rnTrn++ if ( $A[1] eq "SEQ" );

        if ( $ins > 0 ) {
            $rnIns++;
            $rSumIns += $ins;
        }
    }
    FileClose($fhi, $OPT_DiffRFile);
    
    #-- Count query diff features and indels
    $fhi = FileOpen("<", $OPT_DiffQFile);
    while (<$fhi>) {
        chomp;
        my @A = split "\t";
        defined($A[4])
            or die "ERROR: Unrecognized format $OPT_DiffRFile, aborting.\n";
        my $gap = $A[4];
        my $ins = $gap;

        #-- Check for tandem signature
        if ( $A[1] eq "GAP" ) {
            scalar(@A) == 7
                or die "ERROR: Unrecognized format $OPT_DiffRFile, aborting.\n";
            $ins = $A[6] if ( $A[6] > $gap );
            if ( $A[4] <= 0 && $A[5] <= 0 && $A[6] > 0 ) {
                $qnTIns++;
                $qSumTIns += $A[6];
            }
        }

        $qnABases -= $gap if ( $gap > 0 );
        $qnInv++ if ( $A[1] eq "INV" );
        $qnRel++ if ( $A[1] eq "JMP" );
        $qnTrn++ if ( $A[1] eq "SEQ" );

        if ( $ins > 0 ) {
            $qnIns++;
            $qSumIns += $ins;
        }
    }
    FileClose($fhi, $OPT_DiffQFile);

    #-- Count SNPs
    $fhi = FileOpen("<", $OPT_SnpsFile);
    while(<$fhi>) {
        chomp;
        my @A = split "\t";
        scalar(@A) == 12
            or die "ERROR: Unrecognized format $OPT_SnpsFile, aborting\n";

        my $r = uc($A[1]);
        my $q = uc($A[2]);

        $rqSNPs{$r}{$q}++;
        if ( $r eq '.' || $q eq '.' ) { $rqnIndels++; }
        else                          { $rqnSNPs++; }

        if ( $A[4] >= $SNPBuff ) {
            $rqGSNPs{$r}{$q}++;
            if ( $r eq '.' || $q eq '.' ) { $rqnGIndels++; }
            else                          { $rqnGSNPs++; }
        }
    }
    FileClose($fhi, $OPT_SnpsFile);


    #-- Output report
    $fho = FileOpen(">", $OPT_ReportFile);

    print  $fho $header;
    printf $fho "%-15s %20s %20s\n", "", "[REF]", "[QRY]";

    print  $fho "[Sequences]\n";

    printf $fho "%-15s %20d %20d\n",
    "TotalSeqs", $rnSeqs, $qnSeqs;
    printf $fho "%-15s %20s %20s\n",
    "AlignedSeqs",
    ( sprintf "%10d(%.2f%%)",
      $rnASeqs, ($rnSeqs ? $rnASeqs / $rnSeqs * 100.0 : 0) ),
    ( sprintf "%10d(%.2f%%)",
      $qnASeqs, ($rnSeqs ? $qnASeqs / $qnSeqs * 100.0 : 0) );
    printf $fho "%-15s %20s %20s\n",
    "UnalignedSeqs",
     ( sprintf "%10d(%.2f%%)",
       $rnSeqs - $rnASeqs,
       ($rnSeqs ? ($rnSeqs - $rnASeqs) / $rnSeqs * 100.0 : 0) ),
     ( sprintf "%10d(%.2f%%)",
       $qnSeqs - $qnASeqs,
       ($qnSeqs ? ($qnSeqs - $qnASeqs) / $qnSeqs * 100.0 : 0) );
     
    print  $fho "\n[Bases]\n";

    printf $fho "%-15s %20d %20d\n",
    "TotalBases", $rnBases, $qnBases;
    printf $fho "%-15s %20s %20s\n",
    "AlignedBases",
    ( sprintf "%10d(%.2f%%)",
      $rnABases, ($rnBases ? $rnABases / $rnBases * 100.0 : 0) ),
    ( sprintf "%10d(%.2f%%)",
      $qnABases, ($qnBases ? $qnABases / $qnBases * 100.0 : 0) );
    printf $fho "%-15s %20s %20s\n",
    "UnalignedBases",
    ( sprintf "%10d(%.2f%%)",
      $rnBases - $rnABases,
      ($rnBases ? ($rnBases - $rnABases) / $rnBases * 100.0 : 0) ),
    ( sprintf "%10d(%.2f%%)",
      $qnBases - $qnABases,
      ($qnBases ? ($qnBases - $qnABases) / $qnBases * 100.0 : 0) );

    print  $fho "\n[Alignments]\n";

    printf $fho "%-15s %20d %20d\n",
    "1-to-1", $rqnAlignsO, $rqnAlignsO;
    printf $fho "%-15s %20d %20d\n",
    "TotalLength", $rSumLenO, $qSumLenO;
    printf $fho "%-15s %20.2f %20.2f\n",
    "AvgLength",
    ($rqnAlignsO ? $rSumLenO / $rqnAlignsO : 0),
    ($rqnAlignsO ? $qSumLenO / $rqnAlignsO : 0);
    printf $fho "%-15s %20.2f %20.2f\n",
    "AvgIdentity",
    ($rqnAlignsO ? $rqSumIdyO / $rqnAlignsO : 0),
    ($rqnAlignsO ? $rqSumIdyO / $rqnAlignsO : 0);

    print  $fho "\n";

    printf $fho "%-15s %20d %20d\n",
    "M-to-M", $rqnAlignsM, $rqnAlignsM;
    printf $fho "%-15s %20d %20d\n",
    "TotalLength", $rSumLenM, $qSumLenM;
    printf $fho "%-15s %20.2f %20.2f\n",
    "AvgLength",
    ($rqnAlignsM ? $rSumLenM / $rqnAlignsM : 0),
    ($rqnAlignsM ? $qSumLenM / $rqnAlignsM : 0);
    printf $fho "%-15s %20.2f %20.2f\n",
    "AvgIdentity",
    ($rqnAlignsM ? $rqSumIdyM / $rqnAlignsM : 0),
    ($rqnAlignsM ? $rqSumIdyM / $rqnAlignsM : 0);

    print  $fho "\n[Feature Estimates]\n";

    printf $fho "%-15s %20d %20d\n",
    "Breakpoints", $rnBrk, $qnBrk;
    printf $fho "%-15s %20d %20d\n",
    "Relocations", $rnRel, $qnRel;
    printf $fho "%-15s %20d %20d\n",
    "Translocations", $rnTrn, $qnTrn;
    printf $fho "%-15s %20d %20d\n",
    "Inversions", $rnInv, $qnInv;

    print  $fho "\n";

    printf $fho "%-15s %20d %20d\n",
    "Insertions", $rnIns, $qnIns;
    printf $fho "%-15s %20d %20d\n",
    "InsertionSum", $rSumIns, $qSumIns;
    printf $fho "%-15s %20.2f %20.2f\n",
    "InsertionAvg",
    ($rnIns ? $rSumIns / $rnIns : 0),
    ($qnIns ? $qSumIns / $qnIns : 0);

    print  $fho "\n";

    printf $fho "%-15s %20d %20d\n",
    "TandemIns", $rnTIns, $qnTIns;
    printf $fho "%-15s %20d %20d\n",
    "TandemInsSum", $rSumTIns, $qSumTIns;
    printf $fho "%-15s %20.2f %20.2f\n",
    "TandemInsAvg",
    ($rnTIns ? $rSumTIns / $rnTIns : 0),
    ($qnTIns ? $qSumTIns / $qnTIns : 0);

    print  $fho "\n[SNPs]\n";

    printf $fho "%-15s %20d %20d\n",
    "TotalSNPs", $rqnSNPs, $rqnSNPs;
    foreach my $r ("A","C") {
        foreach my $q ("A","C","G","T") {
            if ( $r ne $q ) {
                printf $fho "%-15s %20s %20s\n",
                "$r$q",
                ( sprintf "%10d(%.2f%%)",
                  $rqSNPs{$r}{$q},
                  ($rqnSNPs ? $rqSNPs{$r}{$q} / $rqnSNPs * 100.0 : 0) ),
                ( sprintf "%10d(%.2f%%)",
                  $rqSNPs{$q}{$r},
                  ($rqnSNPs ? $rqSNPs{$q}{$r} / $rqnSNPs * 100.0 : 0) );
            }
        }
    }

    print  $fho "\n";

    printf $fho "%-15s %20d %20d\n",
    "TotalGSNPs", $rqnGSNPs, $rqnGSNPs;
    foreach my $r ("A","C") {
        foreach my $q ("A","C","G","T") {
            if ( $r ne $q ) {
                printf $fho "%-15s %20s %20s\n",
                "$r$q",
                ( sprintf "%10d(%.2f%%)",
                  $rqGSNPs{$r}{$q},
                  ($rqnGSNPs ? $rqGSNPs{$r}{$q} / $rqnGSNPs * 100.0 : 0) ),
                ( sprintf "%10d(%.2f%%)",
                  $rqGSNPs{$q}{$r},
                  ($rqnGSNPs ? $rqGSNPs{$q}{$r} / $rqnGSNPs * 100.0 : 0) );
            }
        }
    }

    print  $fho "\n";

    printf $fho "%-15s %20d %20d\n",
    "TotalIndels", $rqnIndels, $rqnIndels;
    foreach my $r ("A","C","G","T") {
        foreach my $q (".") {
            if ( $r ne $q ) {
                printf $fho "%-15s %20s %20s\n",
                "$r$q",
                ( sprintf "%10d(%.2f%%)",
                  $rqSNPs{$r}{$q},
                  ($rqnIndels ? $rqSNPs{$r}{$q} / $rqnIndels * 100.0 : 0) ),
                ( sprintf "%10d(%.2f%%)",
                  $rqSNPs{$q}{$r},
                  ($rqnIndels ? $rqSNPs{$q}{$r} / $rqnIndels * 100.0 : 0) );
            }
        }
    }

    print  $fho "\n";

    printf $fho "%-15s %20d %20d\n",
    "TotalGIndels", $rqnGIndels, $rqnGIndels;
    foreach my $r ("A","C","G","T") {
        foreach my $q (".") {
            if ( $r ne $q ) {
                printf $fho "%-15s %20s %20s\n",
                "$r$q",
                ( sprintf "%10d(%.2f%%)",
                  $rqGSNPs{$r}{$q},
                  ($rqnGIndels ? $rqGSNPs{$r}{$q} / $rqnGIndels * 100.0 : 0) ),
                ( sprintf "%10d(%.2f%%)",
                  $rqGSNPs{$q}{$r},
                  ($rqnGIndels ? $rqGSNPs{$q}{$r} / $rqnGIndels * 100.0 : 0) );
            }
        }
    }

    FileClose($fho, $OPT_ReportFile);
}


#--------------------------------------------------------------- FastaSizes ----
# Compute lengths for a multi-fasta file and store in hash reference
sub FastaSizes($$)
{

    my $file = shift;
    my $href = shift;
    my ($tag, $len);

    my $fhi = FileOpen("<", $file);
    while (<$fhi>) {
        chomp;

        if ( /^>/ ) {
            $href->{$tag} = $len if defined($tag);
            ($tag) = /^>(\S+)/;
            $len = 0;
        } else {
            if ( /\s/ ) {
                die "ERROR: Whitespace found in FastA $file, aborting.\n";
            }
            $len += length;
        }
    }
    $href->{$tag} = $len if defined($tag);
    FileClose($fhi, $file);
}


#----------------------------------------------------------------- FileOpen ----
# Open file, return filehandle, or die
sub FileOpen($$)
{
    my ($mode, $name) = @_;
    my $fhi;
    open($fhi, $mode, $name)
        or die "ERROR: Could not open $name, aborting. $!\n";
    return $fhi;
}


#---------------------------------------------------------------- FileClose ----
# Close file, or die
sub FileClose($$)
{
    my ($fho, $name) = @_;
    close($fho) or die "ERROR: Could not close $name, aborting. $!\n"
}


#------------------------------------------------------------------- GetOpt ----
# Get command options and check file permissions
sub GetOpt()
{
    #-- Initialize TIGR::Foundation
    $TIGR = new TIGR::Foundation;
    if ( !defined($TIGR) ) {
        print STDERR "ERROR: TIGR::Foundation could not be initialized";
        exit(1);
    }
    
    #-- Set help and usage information
    $TIGR->setHelpInfo($HELP_INFO);
    $TIGR->setUsageInfo($USAGE_INFO);
    $TIGR->setVersionInfo($VERSION_INFO);
    $TIGR->addDependInfo(@DEPEND_INFO);

    #-- Get options
    my $err = !$TIGR->TIGR_GetOptions
        (
         "d|delta=s"  => \$OPT_DeltaFile,
         "p|prefix=s" => \$OPT_Prefix,
         );
    
    #-- Check if the parsing was successful
    if ( $err
         || (defined($OPT_DeltaFile) && scalar(@ARGV) != 0)
         || (!defined($OPT_DeltaFile) && scalar(@ARGV) != 2) ) {
        $TIGR->printUsageInfo();
        print STDERR "Try '$0 -h' for more information.\n";
        exit(1);
    }

    my @errs;

    $TIGR->isExecutableFile($DELTA_FILTER)
        or push(@errs, $DELTA_FILTER);

    $TIGR->isExecutableFile($SHOW_DIFF)
        or push(@errs, $SHOW_DIFF);

    $TIGR->isExecutableFile($SHOW_SNPS)
        or push(@errs, $SHOW_SNPS);

    $TIGR->isExecutableFile($SHOW_COORDS)
        or push(@errs, $SHOW_COORDS);

    $TIGR->isExecutableFile($NUCMER)
        or push(@errs, $NUCMER);

    if ( defined($OPT_DeltaFile) ) {
        $TIGR->isReadableFile($OPT_DeltaFile)
            or push(@errs, $OPT_DeltaFile);

        my $fhi = FileOpen("<", $OPT_DeltaFile);
        $_ = <$fhi>;
        FileClose($fhi, $OPT_DeltaFile);

        ($OPT_RefFile, $OPT_QryFile) = /^(.+) (.+)$/;
    }
    else {
        $OPT_RefFile = File::Spec->rel2abs($ARGV[0]);
        $OPT_QryFile = File::Spec->rel2abs($ARGV[1]);
    }

    $TIGR->isReadableFile($OPT_RefFile)
        or push(@errs, $OPT_RefFile);

    $TIGR->isReadableFile($OPT_QryFile)
        or push(@errs, $OPT_QryFile);

    $OPT_DeltaFileO = $OPT_Prefix . $OPT_DeltaFileO;
    $TIGR->isCreatableFile("$OPT_DeltaFileO")
        or $TIGR->isWritableFile("$OPT_DeltaFileO")
        or push(@errs, "$OPT_DeltaFileO");

    $OPT_DeltaFileM = $OPT_Prefix . $OPT_DeltaFileM;
    $TIGR->isCreatableFile("$OPT_DeltaFileM")
        or $TIGR->isWritableFile("$OPT_DeltaFileM")
        or push(@errs, "$OPT_DeltaFileM");

    $OPT_CoordsFileO = $OPT_Prefix . $OPT_CoordsFileO;
    $TIGR->isCreatableFile("$OPT_CoordsFileO")
        or $TIGR->isWritableFile("$OPT_CoordsFileO")
        or push(@errs, "$OPT_CoordsFileO");

    $OPT_CoordsFileM = $OPT_Prefix . $OPT_CoordsFileM;
    $TIGR->isCreatableFile("$OPT_CoordsFileM")
        or $TIGR->isWritableFile("$OPT_CoordsFileM")
        or push(@errs, "$OPT_CoordsFileM");

    $OPT_SnpsFile = $OPT_Prefix . $OPT_SnpsFile;
    $TIGR->isCreatableFile("$OPT_SnpsFile")
        or $TIGR->isWritableFile("$OPT_SnpsFile")
        or push(@errs, "$OPT_SnpsFile");

    $OPT_DiffRFile = $OPT_Prefix . $OPT_DiffRFile;
    $TIGR->isCreatableFile("$OPT_DiffRFile")
        or $TIGR->isWritableFile("$OPT_DiffRFile")
        or push(@errs, "$OPT_DiffRFile");

    $OPT_DiffQFile = $OPT_Prefix . $OPT_DiffQFile;
    $TIGR->isCreatableFile("$OPT_DiffQFile")
        or $TIGR->isWritableFile("$OPT_DiffQFile")
        or push(@errs, "$OPT_DiffQFile");

    $OPT_ReportFile = $OPT_Prefix . $OPT_ReportFile;
    $TIGR->isCreatableFile("$OPT_ReportFile")
        or $TIGR->isWritableFile("$OPT_ReportFile")
        or push(@errs, "$OPT_ReportFile");

    if ( scalar(@errs) ) {
        print STDERR "ERROR: The following critical files could not be used\n";
        while ( scalar(@errs) ) { print(STDERR pop(@errs),"\n"); }
        print STDERR "Check your paths and file permissions and try again\n";
        exit(1);
    }
}
