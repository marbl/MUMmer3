#!__PERL_PATH

#-------------------------------------------------------------------------------
#   Programmer: Adam M Phillippy, The Institute for Genomic Research
#         File: mummerplot
#         Date: 01 / 08 / 03
#
#        Usage:
#    mummerplot  [options]  <match file>
#
#                Try 'mummerplot -h' for more information.
#
#      Purpose: To create a gnuplot script and datafile to generate a dotplot
#               of the alignments created by mummer, nucmer or promer.
#
#-------------------------------------------------------------------------------

use lib "__SCRIPT_DIR";
use Foundation;
use strict;

my $SCRIPT_DIR = "__SCRIPT_DIR";

my @DEPEND_INFO =
    (
     "$SCRIPT_DIR/Foundation.pm"
     );

my $HELP_INFO = q~
  USAGE: mummerplot  [options]  <match file>

  DESCRIPTION:
    mummerplot creates a gnuplot script and a datafile to generate a dotplot
    of the alignment data contained in the output of mummer, nucmer or promer.
    These files will be output in out.gp, and out.plot. Simply type
    'gnuplot out.gp' to view the plot, or edit the out.gp file to change
    output file, line thickness, labels, color, etc. Explore the out.plot
    file to see the data point collection. *Note all inputs except for a list
    of MUMs will be scaled by sequence length, a plot of MUMs will be auto-
    scaled based on the data points.

    mummerplot will automatically identify the format of the <match file> it is
    given. It accepts a (3 column ONLY) mum list from mummer, a cluster file
    from nucmer or promer, a delta file from nucmer or promer, or a tiling from
    show-tiling.

  MANDATORY:
    match file      Set the alignment input file to "match file"
                    Valid types are (.out, .cluster, .delta and .tiling)

  OPTIONS:
    -c
    --coverage      Generate a reference coverage plot (default for tiling)
    --[no]color     Toggle color option (does not effect X11 plots)
    --depend        Print the dependency information and exit
    -h
    --help          Display help information and exit
    -p|prefix       Set the prefix of the output files (default 'out')
    --postscript    Set the plot terminal to postscript
    -r|IdR          Use reference sequence ID for the X-axis sequence
                    Useful for selecting one reference sequence for plotting,
                    has no effect for mummer2 input
    -q|IdQ          Use query sequence ID for the Y-axis sequence
                    Useful for selecting one query sequence for plotting
    -x|xrange       Set the xrange for the plot "[min,max]"
    -y|yrange       Set the yrange for the plot "[min,max]"
    -V
    --version       Display the version information and exit
    ~;


my $USAGE_INFO = q~
  USAGE: mummerplot  [options]  <match file>
    ~;


my $VERSION_INFO = q~
mummerplot version 2.0
    ~;


my $MUMMER_STRING = "MUMMER";
my $CLUSTER_STRING = "CLUSTER";
my $DELTA_STRING = "DELTA";
my $TILING_STRING = "TILING";

my $X11_TERM = "X11";
my $POSTSCRIPT_TERM = "postscript";


#-- Global for efficiency, not elegance
my $StartRef;         # Start of match in reference
my $StartQry;         # Start of match in query
my $Sim;              # Approx %similarity
my $LenRef;           # Length of match in reference
my $LenQry;           # Length of match in query
my $DirRef;           # Direction of the reference
my $DirQry;           # Direction of the query
my $DataWasPlotted;   # Was data plotted?
my $CoveragePlot;     # Is this a one vs many plot?
my $xmin;             # Plot ranges
my $xmax;
my $ymin;
my $ymax;


sub main ( )
{
    my $tigr;             # TIGR::Foundation object
    my $err;              # Error variable

    my $matchfilename;    # Input file name
    my $scriptfilename;   # gnuplot script file
    my $fwddatafilename;  # gnuplot forward match data file
    my $revdatafilename;  # gnuplot reverse match data file
    my $datafilename;     # single data file to replace the above two

    my $format;           # Input file format
    my $prefix = "out";   # Output file prefix
    my $color = 1;        # color plot?

    my $multRef;          # multiple references were plotted
    my $multQry;          # multiple queries were plotted

    my $currlenRefSeq;        # Length of the refernece sequence
    my $currlenQrySeq;        # Length of the query sequence
    my $lenRefSeq;
    my $lenQrySeq;

    my $prevRefID;        # Current reference sequence ID
    my $prevQryID;        # Current query sequence ID
    my $currRefID;        # Current reference sequence ID
    my $currQryID;        # Current query sequence ID
    my $useRefID;         # Use reference sequence ID
    my $useQryID;         # Use query sequence ID

    my $foundRefID;       # was the reference ID found?
    my $foundQryID;       # was the query ID found?

    my $endRef;           # End of match in reference
    my $endQry;           # End of match in query

    my $xrange;
    my $yrange;

    my $terminal;

    my $line;
    my $delta;
    my $total;
    my $datatype;

    #-- Initialize TIGR::Foundation
    $tigr = new TIGR::Foundation;
    if ( !defined ($tigr) ) {
        print (STDERR "ERROR: TIGR::Foundation could not be initialized");
        exit (1);
    }
  
    #-- Set help and usage information
    $tigr->setHelpInfo ($HELP_INFO);
    $tigr->setUsageInfo ($USAGE_INFO);
    $tigr->setVersionInfo ($VERSION_INFO);
    $tigr->addDependInfo (@DEPEND_INFO);

    #-- Get command line parameters
    $err = $tigr->TIGR_GetOptions
        (
	 "c|coverage" => \$CoveragePlot,
	 "color!" => \$color,
	 "p|prefix=s" => \$prefix,
	 "r|IdR=s" => \$useRefID,
	 "q|IdQ=s" => \$useQryID,
	 "postscript" => \$terminal,
	 "x|xrange=s" => \$xrange,
	 "y|yrange=s" => \$yrange
         );

    #-- Check if the parsing was successful
    if ( $err == 0  ||  scalar(@ARGV) != 1 ) {
        $tigr->printUsageInfo( );
        print (STDERR "Try '$0 -h' for more information.\n");
        exit (1);
    }

    if ( defined ($xrange) ) {
	if ( ! (($xmin,$xmax) = $xrange =~ /\[(\d+),(\d+)\]/) ) {
	    die "ERROR: invalid xrange\n";
	}
    }

    if ( defined ($yrange) ) {
	if ( ! (($ymin,$ymax) = $yrange =~ /\[(\d+),(\d+)\]/) ) {
	    die "ERROR: invalid yrange\n";
	}
    }

    if ( defined ($terminal) ) {
	$terminal = $POSTSCRIPT_TERM;
    } else {
	$terminal = $X11_TERM;
    }

    #-- Open input file
    $matchfilename = $ARGV[0];
    open(MATCH_FILE, "<$matchfilename")
	or die "ERROR: could not open $matchfilename $!\n";


    #-- Open output files
    $scriptfilename = $prefix . ".gp";
    open(SCRIPT_FILE, ">$scriptfilename")
	or die "ERROR: could not open $scriptfilename $!\n";
    $fwddatafilename = $prefix . ".fwdplot." . $$;
    open(FWDDATA_FILE, "+>$fwddatafilename")
	or die "ERROR: could not open $fwddatafilename $!\n";
    $revdatafilename = $prefix . ".revplot." . $$;
    open(REVDATA_FILE, "+>$revdatafilename")
	or die "ERROR: could not open $revdatafilename $!\n";
    $datafilename = $prefix . ".plot";
    open(DATA_FILE, ">$datafilename")
	or die "ERROR: could not open $datafilename $!\n";


    #-- Determine input file format
    chomp ( $line = <MATCH_FILE> );
    if ( $line =~ /^>\S+\s+\d+\s+bases/ ) {
	#-- It'a a tiling output file
	
	$format = $TILING_STRING;
	$CoveragePlot = 1;
	$currQryID = "";
	($currRefID,$currlenRefSeq) = $line =~ /^>(\S+)\s+(\d+)/;

    } elsif ( $line =~ /^>\s+\S+/ ) {
	#-- It's a mummer output file

	$format = $MUMMER_STRING;

	$useRefID = undef;
	$currRefID = "";
	($currQryID) = $line =~ /^>\s*(\S+)/;
	if ( $line =~ /^>\s*\S+\s+Reverse$/ ) {
	    $DirRef = 1;
	    $DirQry = -1;
	} else {
	    $DirRef = 1;
	    $DirQry = 1;
	}
    } else {
	#-- It's a nucmer/promer output file

	chomp ($line = <MATCH_FILE>);
	if ( !($line =~ /^NUCMER$/)  &&
	     !($line =~ /^PROMER$/) ) {
	    die "ERROR: $matchfilename is unrecognized file type\n" .
		"valid file types are (.out, .cluster, .delta, .tiling)\n" .
		    "suspicious line: $line\n";
	}
	$datatype = $line;
	chomp ($line = <MATCH_FILE>);
	if ( $line =~ /^>\S+\s+\S+\s+\d+\s+\d+$/ ) {
	    $format = undef( );
	} else {
	    die "ERROR: $matchfilename is unrecognized file type\n" .
		"valid file types are (.out, .cluster, .delta, .tiling)\n" .
		    "suspicious line: $line\n";
	}

	($currRefID, $currQryID, $currlenRefSeq, $currlenQrySeq) =
	    $line =~ /^>(\S+)\s+(\S+)\s+(\d+)\s+(\d+)/;
    }


    if ( defined ($format) ) {
	print STDERR "INFO: assuming $format format\n";
	if ( $format eq $MUMMER_STRING ) {
	    print STDERR
	     "INFO: the '-b' and '-c' options must have been used when\n" .
	     "      running mummer, or the resulting plot may be incorrect\n";
	}
    }


    while ( $line = <MATCH_FILE> ) {
	chomp $line;

	#-- Make sure the sequence names are defined
	if ( !defined($currRefID)  ||  !defined($currQryID) ) {
	    die "ERROR: $matchfilename is unrecognized file type\n";
	}


	if ( $line =~ /^>/ ) {
	    #-- It's a sequence header line

	    $prevRefID = $currRefID;
	    $prevQryID = $currQryID;

	    #-- Read the new headers
	    if ( $format eq $TILING_STRING ) {
		($currRefID,$currlenRefSeq) = $line =~ /^>(\S+)\s+(\d+)/;

	    } elsif ( $format eq $MUMMER_STRING ) {
		($currQryID) = $line =~ /^>\s*(\S+)/;

		#-- Set the direction of the sequences
		if ( $line =~ /^>\s*\S+\s+Reverse$/ ) {
		    $DirRef = 1;
		    $DirQry = -1;
		} else {
		    $DirRef = 1;
		    $DirQry = 1;
		}
	    } else {
		($currRefID, $currQryID, $currlenRefSeq, $currlenQrySeq) =
		    $line =~ /^>(\S+)\s+(\S+)\s+(\d+)\s+(\d+)/;
	    }


	    #-- Warn if multiple sequences are used in the same plot
	    if ( !$multRef  &&  $prevRefID ne $currRefID  &&
		 !defined ($useRefID) ) {
		print STDERR
		  "WARNING: multiple reference sequences in the same plot,\n" .
		  "         use the '-r ID' option to aviod confusion\n";
		$multRef = 1;
	    }
	    if ( !$multQry  &&  $prevQryID ne $currQryID  &&
		 !defined ($useQryID) ) {
		if ( !$CoveragePlot ) {
		  print STDERR
		    "WARNING: multiple query sequences in the same plot,\n" .
		    "         use the '-q ID' option to aviod confusion\n";
		}
		$multQry = 1;
	    }

	} else {
	    #-- It's a data line

	    #-- If not the selected sequence, continue
	    if ( defined ($useRefID)  &&  $useRefID ne $currRefID ) {
		if ( defined ($useQryID)  &&  $useQryID ne $currQryID ) {
		    next;
		}
		$foundQryID = 1;
		next;
	    } else {
		$foundRefID = 1;
		if ( defined ($useQryID)  &&  $useQryID ne $currQryID ) {
		    next;
		}
		$foundQryID = 1;
	    }

	    $lenRefSeq = $currlenRefSeq;
	    $lenQrySeq = $currlenQrySeq;


	    #-- If you haven't figured it out yet, check type
	    if ( !defined ($format) ) {
		if ( $line =~ /^[\s\-][1-3] [\s\-][1-3]$/ ) {
		    $format = $CLUSTER_STRING;
		} elsif ( $line =~
			  /^\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+$/ ) {
		    $format = $DELTA_STRING;
		} else {
		    die "ERROR: $matchfilename is unrecognized file type\n" .
		 "valid file types are (.out, .cluster, .delta, .tiling)\n" .
		 "suspicious line: $line\n";
		}
		print STDERR "INFO: assuming $format format\n";
	    }


	    if ( $format eq $TILING_STRING ) {
		#-- It's a tiling match line
		($StartRef,$_,$_,$LenRef,$_,$Sim,$DirQry,$currQryID) =
		    split " ", $line;
		$DirRef = 1;
		$StartQry = 0;
		$LenQry = $LenRef;
		if ( $DirQry eq "+" ) {
		    $DirQry = 1;
		} elsif ( $DirQry eq "-" ) {
		    $DirQry = -1;
		} else {
		    die "ERROR: Invalid tile direction \"$DirQry\"\n";
		}

		#-- Output the tile to the plot files
		outputCurrentMatch ( );

	    } elsif ( $format eq $MUMMER_STRING ) {
		#-- It's a mummer match line
		
		if ( ! ($line =~ /^\s*\d+\s*\d+\s*\d+$/) ) {
		    die "ERROR: $matchfilename is unrecognized file type\n" .
		 "valid file types are (.out, .cluster, .delta, .tiling)\n" .
		 "suspicious line: $line\n";
		}
		
		#-- Pull the data
		($StartRef,$StartQry,$LenRef) = split " ", $line;
		$LenQry = $LenRef;
		
		#-- Output the match to the plot files
		outputCurrentMatch ( );

	    } elsif ( $format eq $CLUSTER_STRING ) {
		#-- It's a cluster header or match line
		
		if ( $line =~ /^[\s\-][1-3] [\s\-][1-3]$/ ) {
		    #-- It's a cluster header
		    
		    #-- Set the direction of the sequences
		    if ( $line =~ /^-/ ) {
			$DirRef = -1;
		    } else {
			$DirRef = 1;
		    }
		    if ( $line =~ /-[1-3]$/ ) {
			$DirQry = -1;
		    } else {
			$DirQry = 1;
		    }
		} elsif ( $line =~
			  /^\s*\d+\s+\d+\s+\d+\s+[\d\-]+\s+[\d\-]+$/ ) {
		    #-- It's a cluster data line
		    
		    #-- Pull the data
		    ($StartRef,$StartQry,$LenRef) = split " ", $line;
		    $LenQry = $LenRef;
		    
		    #-- Output the match to the plot files
		    outputCurrentMatch ( );
		} else {
		    die "ERROR: $matchfilename is unrecognized file type\n" .
	      	 "valid file types are (.out, .cluster, .delta, .tiling)\n" .
	    	 "suspicious line: $line\n";
		}
		
	    } elsif ( $format eq $DELTA_STRING ) {
		#-- It's an alignment header or delta-int
		
		if ( $line =~ /^[-]?\d+$/ ) {
		    #-- It's a delta-int
		    
		    ($delta) = $line =~ /(\S+)/;
		    if ( $delta < 0 ) {
			$total ++;
		    } elsif ( $delta == 0 ) {
			$Sim = ($total - $Sim) / $total * 100.0;

			#-- Output the match to the plot files
			outputCurrentMatch ( );
		    }
		    
		} elsif ( $line =~
			  /^\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+$/ ) {
		    #-- It's an alignment header
		    
		    #-- Pull the data
		    ($StartRef,$endRef,$StartQry,$endQry,$_,$Sim) =
			split " ", $line;
		    $LenRef = abs($endRef - $StartRef) + 1;
		    $LenQry = abs($endQry - $StartQry) + 1;

		    $total = $LenRef;
		    if ( $datatype =~ /^PROMER$/ ) {
			$total /= 3.0;
		    }

		    #-- Set the direction of the sequences
		    if ( $StartRef > $endRef ) {
			$DirRef = -1;
		    } else {
			$DirRef = 1;
		    }
		    if ( $StartQry > $endQry ) {
			$DirQry = -1;
		    } else {
			$DirQry = 1;
		    }

		} else {
		    die "ERROR: $matchfilename is unrecognized file type\n" .
		 "valid file types are (.out, .cluster, .delta, .tiling)\n" .
		 "suspicious line: $line\n";
		}
		
	    } else {
		die "ERROR: invalid format type '$format'\n";
	    }
	}
    }


    close(MATCH_FILE)
	or die "ERROR: could not close $matchfilename $!\n";


    print DATA_FILE "#-- FORWARD MATCHES\n";
    print DATA_FILE "0 0\n\n";
    seek FWDDATA_FILE, 0, 0;
    while ( <FWDDATA_FILE> ) {
	print DATA_FILE;
    }
    print DATA_FILE "\n#-- REVERSE MATCHES\n";
    print DATA_FILE "0 0\n\n";
    seek REVDATA_FILE, 0, 0;
    while ( <REVDATA_FILE> ) {
	print DATA_FILE;
    }

    close(FWDDATA_FILE)
	or die "ERROR: could not close $fwddatafilename $!\n";
    close(REVDATA_FILE)
	or die "ERROR: could not close $revdatafilename $!\n";
    close(DATA_FILE)
	or die "ERROR: could not close $datafilename $!\n";


    $err = unlink $revdatafilename, $fwddatafilename;
    if ( $err != 2 ) {
	print STDERR "WARNING: could not delete temporary files,\n" .
	    "$revdatafilename and $fwddatafilename\n";
    }


    #-- Dynamically generate the gnuplot script
    printf SCRIPT_FILE "set terminal $terminal %s\n",
         $terminal eq $POSTSCRIPT_TERM ? $color ? "color solid" : "solid" : "";
    if ( $terminal eq $POSTSCRIPT_TERM ) {
	print SCRIPT_FILE "set output \"$prefix.ps\"\n";
    }
    if ( $CoveragePlot  &&  $terminal eq $POSTSCRIPT_TERM ) {
	print SCRIPT_FILE "set size 1,.5\n";
    }
    print SCRIPT_FILE "set nokey\n";
    print SCRIPT_FILE "set nogrid\n";
    printf SCRIPT_FILE "set xlabel %s\n", defined ($useRefID) ?
	"\"$useRefID\"" : $multRef ? "\"References\"" : "\"$currRefID\"";
    printf SCRIPT_FILE "set ylabel %s\n", $CoveragePlot ? "\"%similarity\"" :
	defined ($useQryID) ? "\"$useQryID\"" : $multQry ?	
	"\"Queries\"" : "\"$currQryID\"";
    printf SCRIPT_FILE "set xrange [%s:%s]\n", defined ($xmin) ? $xmin : 0,
    defined ($xmax) ? $xmax :
    defined ($lenRefSeq) && !$multRef ? "$lenRefSeq" : "*";
    printf SCRIPT_FILE "set yrange [%s:%s]\n", defined ($ymin) ? $ymin : 0,
    defined ($ymax) ? $ymax : $CoveragePlot ? "110" :
	defined($lenQrySeq)  &&  !$multQry ? "$lenQrySeq" : "*";
    print SCRIPT_FILE "plot ";
    if ( $CoveragePlot ) {
	print SCRIPT_FILE "\"$datafilename\" index 0 title \"Fwd\" w l lt 1 lw 5, ";
	print SCRIPT_FILE "\"$datafilename\" index 1 title \"Rev\" w l lt 2 lw 5\n";
    } else {
	print SCRIPT_FILE "\"$datafilename\" index 0 title \"Fwd\" w lp lt 1 lw 2 pt 1 ps .75, ";
	print SCRIPT_FILE "\"$datafilename\" index 1 title \"Rev\" w lp lt 2 lw 2 pt 1 ps .75\n";
    }
    if ( $terminal eq $X11_TERM ) {
	print SCRIPT_FILE "pause -1 \"Hit return to continue\"\n";
    }


    close(SCRIPT_FILE)
	or die "ERROR: could not close $scriptfilename $!\n";


    if ( !defined($DataWasPlotted) ) {
	print STDERR "WARNING: no data was plotted\n";
    }
    if ( defined($useRefID)  &&  !defined($foundRefID) ) {
	print STDERR "WARNING: reference sequence $useRefID was not found\n";
    }
    if ( defined($useQryID)  &&  !defined($foundQryID) ) {
	print STDERR "WARNING: query sequence $useQryID was not found\n";
    }
}


sub outputCurrentMatch ( ) {

    my $tmp;

    $DataWasPlotted = 1;

    if ( !(defined($StartRef) && defined($StartQry) &&
	   defined($LenRef) && defined($LenQry) &&
	   defined($DirRef) && defined($DirQry)) ) {
	die "ERROR: attempted to output undefined values\n" .
	    "       please file a bug report\n";
    }

    if ( $DirRef == $DirQry ) {
	#-- Forward

	if ( $DirRef == -1 ) {
	    $StartRef -= $LenRef - 1;
	    $StartQry -= $LenQry - 1;
	}

	printf FWDDATA_FILE "%s\n", $CoveragePlot ?
	    "$StartRef 10" : "$StartRef $StartQry";
	$tmp = $StartRef;
	$StartRef += $LenRef - 1;
	$StartQry += $LenQry - 1;
	printf FWDDATA_FILE "%s\n\n", $CoveragePlot ?
	    "$StartRef 10" : "$StartRef $StartQry";

	if ( $CoveragePlot  &&  defined ($Sim) ) {
	    print FWDDATA_FILE "$tmp $Sim\n";
	    print FWDDATA_FILE "$StartRef $Sim\n\n";
	}
    } else {
	#-- Reverse complement

	printf REVDATA_FILE "%s\n", $CoveragePlot ?
	    "$StartRef 10" : "$StartRef $StartQry";
	$tmp = $StartRef;
	if ( $DirRef == -1 ) {
	    $StartRef -= $LenRef - 1;
	    $StartQry += $LenQry - 1;
	} else { 
	    $StartRef += $LenRef - 1;
	    $StartQry -= $LenQry - 1;
	}
	printf REVDATA_FILE "%s\n\n", $CoveragePlot ?
	    "$StartRef 10" : "$StartRef $StartQry";

	if ( $CoveragePlot  &&  defined ($Sim) ) {
	    print REVDATA_FILE "$tmp $Sim\n";
	    print REVDATA_FILE "$StartRef $Sim\n\n";
	}
    }
}

exit ( main ( ) );

#-- END OF SCRIPT
