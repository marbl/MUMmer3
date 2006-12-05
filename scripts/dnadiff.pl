#!__PERL_PATH

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

my $OPT_Prefix     = "out";         # prefix for all output files
my $OPT_RefFile;                    # reference file
my $OPT_QryFile;                    # query file
my $OPT_DeltaFile;                  # unfiltered alignment file
my $OPT_FilterFile = ".deltaf";     # filtered alignment file
my $OPT_SnpsFile   = ".snp";        # snps output file
my $OPT_DiffFile   = ".diff";       # difffle file
my $OPT_ReportFile = ".report";     # report file

my $TIGR;  # TIGR Foundation object


sub RunNucmer();
sub RunFilter();
sub RunSNPs();

sub GetOpt();


#--------------------------------------------------------------------- main ----
 main:
{
    GetOpt();

    RunNucmer() unless defined($OPT_DeltaFile);
    RunFilter();
    RunSNPs();

    exit(0);
}


#------------------------------------------------------------ MakeAlignment ----
sub MakeAlignment()
{
    print STDERR ("Building alignments\n");
    my $cmd = "$NUCMER --maxmatch -p $OPT_Prefix $OPT_RefFile $OPT_QryFile";
    my $err = "ERROR: failed to run nucmer, aborting.\n";

    system($cmd) == 0 or die $err;
    $OPT_DeltaFile = $OPT_Prefix . ".delta";
}


#---------------------------------------------------------------- RunFilter ----
sub RunFilter()
{
    print STDERR ("Filtering alignments\n");
    my $cmd = "$DELTA_FILTER -r -q $OPT_DeltaFile > $OPT_FilterFile";
    my $err = "ERROR: failed to run delta-filter, aborting.\n";

    system($cmd) == 0 or die $err;
}


#------------------------------------------------------------------ RunSNPs ----
sub RunSNPs()
{
    print STDERR ("Analyzing SNPs\n");
    my $cmd = "$SHOW_SNPS -rlTC $OPT_FilterFile > $OPT_SnpsFile";
    my $err = "ERROR: failed to run show-snps, aborting.\n";

    system($cmd) == 0 or die $err;
}


#------------------------------------------------------------------- GetOpt ----
sub GetOpt()
{
    #-- Initialize TIGR::Foundation
    $TIGR = new TIGR::Foundation;
    if ( !defined($TIGR) ) {
        print STDERR ("ERROR: TIGR::Foundation could not be initialized");
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
    if ( $err  ||  scalar(@ARGV) != 2 ) {
        $TIGR->printUsageInfo();
        print STDERR ("Try '$0 -h' for more information.\n");
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

    $OPT_RefFile = File::Spec->rel2abs ($ARGV[0]);
    $TIGR->isReadableFile($OPT_RefFile)
        or push(@errs, $OPT_RefFile);
    
    $OPT_QryFile = File::Spec->rel2abs ($ARGV[1]);
    $TIGR->isReadableFile($OPT_QryFile)
        or push(@errs, $OPT_QryFile);

    if ( defined($OPT_DeltaFile) ) {
        $TIGR->isReadableFile($OPT_DeltaFile)
            or push(@errs, $OPT_DeltaFile);
    }

    $OPT_FilterFile = $OPT_Prefix . $OPT_FilterFile;
    $TIGR->isCreatableFile("$OPT_FilterFile")
        or $TIGR->isWritableFile("$OPT_FilterFile")
        or push(@errs, "$OPT_FilterFile");

    $OPT_SnpsFile = $OPT_Prefix . $OPT_SnpsFile;
    $TIGR->isCreatableFile("$OPT_SnpsFile")
        or $TIGR->isWritableFile("$OPT_SnpsFile")
        or push(@errs, "$OPT_SnpsFile");

    $OPT_DiffFile = $OPT_Prefix . $OPT_DiffFile;
    $TIGR->isCreatableFile("$OPT_DiffFile")
        or $TIGR->isWritableFile("$OPT_DiffFile")
        or push(@errs, "$OPT_DiffFile");

    $OPT_ReportFile = $OPT_Prefix . $OPT_ReportFile;
    $TIGR->isCreatableFile("$OPT_ReportFile")
        or $TIGR->isWritableFile("$OPT_ReportFile")
        or push(@errs, "$OPT_ReportFile");

    if ( scalar(@errs) ) {
        print STDERR
            ("ERROR: The following critical files could not be used\n");
        while ( scalar(@errs) ) { print(STDERR pop(@errs),"\n"); }
        print STDERR
            ("Check your paths and file permissions and try again\n");
        exit(1);
    }
}
