use ExtUtils::MakeMaker;
use ExtUtils::MakeMaker::Config;
require ExtUtils::Manifest;
use Term::ANSIColor;
require "perl2wrapper_ergatis.pl";
require "python2wrapper_ergatis.pl";
use strict;

$|++;

## always remake the MANIFEST
ExtUtils::Manifest::mkmanifest();

## let's start by checking for the modules we need
my @modules_to_check = qw( 
    Bio::Seq
    CDB_File
    Class::Struct
    CGI::Session
    Config::IniFiles
    DB_File
    Data::Dumper
    Date::Manip
    File::Basename
    File::Copy
    File::Find
    File::Mirror
    File::Path
    File::Spec
    Getopt::Long
    HTML::Template
    IO::File
    IO::Tee
    IPC::Open3
    LWP::Simple
    Log::Cabin
    Log::Log4perl
    Mail::Mailer
    Math::Combinatorics
    MLDBM
    PerlIO::gzip
    Storable
    URI::Escape
    XML::LibXML
    XML::Parser
    XML::RSS
    XML::Twig
    XML::Writer
);
my @modules_not_found;

print "checking necessary perl modules:\n";

for ( @modules_to_check ) {
    my $formed = $_;
    $formed =~ s|::|/|g;
    $formed =~ s|$|.pm|g;

    print "\t$_\t:";
    eval {require $formed};

    if ( $@ ) {
        print color 'bold red';
        print " not found\n";
        push @modules_not_found, $_;
    } else {
        print color 'bold green';
        print " found\n";
    }

    print color 'reset';
}

if ( scalar @modules_not_found && !($ARGV[0] eq 'force')) {
    print "\ninstall aborted.  the following modules are missing: \n";
    for ( @modules_not_found ) {
        print "\t$_\n";
    }

    exit(1);
}

## make sure that their ExtUtils::MakeMaker is a high enough version
my $makemaker_found = 0;

eval "require ExtUtils::MakeMaker";

if ( ExtUtils::MakeMaker->VERSION < 6.31 ) {
    print color 'bold red';
    print "ExtUtils::MakeMaker >= 6.31 not found\n";
    print color 'reset';
    exit(1);
}

print "\nall modules found, carrying on\n\n";

my $install_base;
my $build_dir = '.';
my $cbuild_dir = 'src';
my $tbuild_dir = 'project_saved_templates';
my $abuild_dir = 'autopipe_package';
my $sdir = 'play_data';
my $software_dir = 'software';

foreach my $arg (@ARGV){
    if ($arg =~ /INSTALL_BASE/) {
		($install_base) = ($arg =~ /INSTALL_BASE=(\S+)/);
    } elsif ( $arg =~ /PREFIX/ ) {
        print "\nERROR: please use INSTALL_BASE instead of PREFIX.  (This requires a current version of ExtUtils::MakeMaker";
        exit(1);
    } elsif ( $arg =~/LIVE_BUILD/){
		$build_dir = "../";
		$cbuild_dir = "../src/c";
		$tbuild_dir = "../templates";
    }
}

my $workflowdocsdir = "$install_base/docs";

if (! -e $workflowdocsdir ) {
    mkdir($workflowdocsdir) || die "failed to create WORKFLOW_DOCS_DIR: $!";
}

my $sampledocsdir = "$install_base/samples";

## workflow and schema docs dirs have always been the same
my $schemadocsdir = $workflowdocsdir;

my($perl_scripts,$perl_wrapperstr) = &get_perl_bins();
my($python_scripts,$python_wrapperstr) = &get_python_bins();
my @script_array = (@$perl_scripts,@$python_scripts);
my $wrapperstr = (join(' ',$perl_wrapperstr,$python_wrapperstr));
my $scripts = \@script_array;

print "Writing makefile to install directory $build_dir to $install_base\n";
my $mm = WriteMakefile(
	      'NAME'=>'ergatis',
	      'EXE_FILES'=>$scripts,
	      macro=>{BUILD_DIR => "$build_dir", A_BUILD_DIR=>"$abuild_dir", T_BUILD_DIR=>"$tbuild_dir", SAMPLE_DIR=>"$sdir", SOFTWARE_DIR=>"$software_dir"},
	      'clean' =>{FILES => "MANIFEST.bak $wrapperstr README.symlinks"}
	      );

package MY;

sub MY::postamble { 
    my $postamble = "WORKFLOW_DOCS_DIR=$workflowdocsdir\nSAMPLE_DOCS_DIR=$sampledocsdir";
    $postamble .= << 'END';

ergatis_install : dirs pipeline_templates autopipe sample_data packages
	@echo "Installing components"
	find $(BUILD_DIR)/components -path 'to_delete' -prune -o -print -type f -not -name software.config -not -name "*.tmpl" -exec install  {} ${WORKFLOW_DOCS_DIR} \;
	#find $(BUILD_DIR)/components -type f -name "*.tmpl" -exec install  {} ${WORKFLOW_DOCS_DIR}/documentation \;
	@install $(BUILD_DIR)/software.config $(PREFIX)

dirs :
	mkdir -p $(WORKFLOW_DOCS_DIR)

pipeline_templates :
	@echo "Installing pipeline templates in $(PREFIX)/project_saved_templates"
	find $(T_BUILD_DIR) -type d -exec mkdir -p $(PREFIX)/{} \;
	find $(T_BUILD_DIR) -path 'to_delete' -prune -o -print -exec install {} $(PREFIX)/{} \;
	find $(T_BUILD_DIR) -path 'to_delete' -prune -o -print -type d -exec install {} ${PREFIX} \;

autopipe :
	@echo "Installing command line pipeline exec in $(PREFIX)/autopipe_package"
	find $(A_BUILD_DIR) -type d -exec mkdir -p $(PREFIX)/{} \;
	find $(A_BUILD_DIR) -path 'to_delete' -prune -o -print -exec install {} $(PREFIX)/{} \;
	find $(A_BUILD_DIR) -path 'to_delete' -prune -o -print -type d -exec install {} ${PREFIX} \;

sample_data :
	@echo "Installing sample data to $(PREFIX)/play_data"
	find $(SAMPLE_DIR) -type d -exec mkdir -p $(PREFIX)/{} \;
	find $(SAMPLE_DIR) -path 'to_delete' -prune -o -print -exec install {} $(PREFIX)/{} \;
	find $(SAMPLE_DIR) -path 'to_delete' -prune -o -print -type d -exec install {} ${PREFIX} \;

packages :
	@echo "Installing pre-packaged software to $(PREFIX)/software"
	find $(SOFTWARE_DIR) -type d -exec mkdir -p $(PREFIX)/{} \;
	find $(SOFTWARE_DIR) -path 'to_delete' -prune -o -print -exec install {} $(PREFIX)/{} \;
	find $(SOFTWARE_DIR) -path 'to_delete' -prune -o -print -type d -exec install {} ${PREFIX} \;


END

    $postamble;
}

sub install { 
    my $self = shift;
    my $installstr = $self->MM::install(@_);
    my $addltarget = 'ergatis_install';
    $installstr =~ s/(pure_install\s+::\s+)(.*)/$1 $2 $addltarget/;
    return $installstr;
}
