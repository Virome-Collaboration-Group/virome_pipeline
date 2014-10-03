#! /usr/local/bin/perl

=head1 NAME

Batch_AAT_Wrap_Move.pl - one-sentence description.

=head1 SYNOPSIS

    USAGE: Batch_AAT_Wrap_Move.pl OPTION LIST HERE

=head1 OPTIONS

=over 4

=item -D / --Database
Annotation database

=item -d / --directory
Repositiory root (default: $ANNOTATION_DIR/$DB/output_repository

=item -p / --pipeline_id
Workflow pipeline ID (it is possible to specify more instances)

=item -L / --id_list
List of pipeline IDs

=item -R / --remove_source
Optional. Remove the results from the source directory after moving them

=item -h / --help
Prints this page and quits

=item -l / --fragment_length
Optional. The length of the fragment searched, if split.  Use 0 for AAT searches
where the input sequence wasn't split.

=back

=head1  DESCRIPTION

I will fill this stuff later

A detailed description of what this program is intended to do, along with current 
limitations.

=head1  INPUT

A description of the input your script requires should go here.  This should include both
a full description of any command-line options as well as input file formats expected.

=head1  OUTPUT

A description of your script's output should go here.

=head1  CONTACT

    Paolo Amedeo
    pamedeo@tigr.org

=begin comment
## legal values for status are active, inactive, hidden, unstable
    status: unstable
    Comment
=end comment

=cut


use strict;
use warnings;
#use lib ('/local/perl/lib/site_perl');
use DBI;
use Getopt::Long;
use Pod::Usage;

MAIN:{
	my ($prN) = ($0 =~ /([^\/]+)$/);
	my $dbfile = "$ENV{EGC_SCRIPTS}/egc_password";
	
	my $safetime = 5; # waiting time before actually running anything where the option -R / remove_source is used
	
	my $debug = 0;
	
	my $mover = '/home/pamedeo/CVS/cram/src/bridge/Move_Workflow_Assembly_Searches.pl'; 
	
	# specific parameters for moving the different kind of results
	
	my %move_opt = (btab => '-B',
	                gap2 => '-A -z',
		        nap  => '-A -z');
                        
        my %tool_tr = (aat_na => 'gap2',
                       aat_aa => 'nap');
	
        my %list_names = (aat_na => {btab => 'aat_na.btab.list', 
                                     gap2 => 'aat_na.gap2.list'},
                          aat_aa => {btab => 'aat_aa.btab.list', 
                                     nap  => 'aat_aa.nap.list'},
                          gap2   => {btab => 'adjust_btab_coordinates.list',
                                     gap2 => 'adjust_gap2_coordinates.list'},
                          nap    => {btab => 'adjust_btab_coordinates.list', 
                                     nap  => 'adjust_nap_coordinates.list'});
	
	my @piplst = ();
	my %pipes = ();
	

	my ($db, $big_db, $id_file, $pipdir, $remove, $help, $fragment_length);
	
	Getopt::Long::Configure("no_ignore_case");
	GetOptions('Database|D=s'        => \$db,
	           'directory|d=s'       => \$pipdir,
		   'pipeline_id|p=i'     => \@piplst,
		   'id_list|L=s'         => \$id_file,
		   'remove_source|R=s'   => \$remove,
		   'fragment_length|l=i' => \$fragment_length,
		   'debug|g=i'           => \$debug,
		   'help|h'              => \$help) || die "\n\nProblems processing the options\n\n";
	

	pod2usage(-exitval => 0,
	          -verbose => 1,
		  -output  => \*STDERR) if $help;
	
	if (defined $db){
		$big_db = uc($db);
	} else {
		pod2usage(-msg     => "\n\nOption -D / --Database (Annotation database) is required\n\n",
		          -exitval => 1,
		          -verbose => 1,
			  -output  => \*STDERR);
	}
	
	if (defined $pipdir && $pipdir =~ /\S/){
		$pipdir =~ s#/$##;
		
		if ($pipdir =~ /$db$/i){
			die "\n\nProblem finding the result directory\n\n" unless -d $pipdir;
		} else {
			if (-d "$pipdir/$big_db"){
				$pipdir .= "/$big_db";
			}
			elsif (-d "$pipdir/$db"){
				$pipdir .= "/$db";
			} else {
				die "\n\nProblem localizing the result directory when looking for $pipdir/$db\n\n";
			}
		}
	} else {
                $pipdir = "$ENV{ANNOTATION_DIR}/$big_db/output_repository";
        }
	
	if (defined $id_file){
		if (open(my $plst, $id_file)){
			while (<$plst>){
				undef $pipes{$1} if /(\d+)/;
			}
			close($plst);
		} else {
			die "\n\nImpossible to open the file $id_file\n\n";
		}
	}
	
	if ($remove){
		print "Warn: you have decided to remove the source files: you have $safetime seconds still to stop this program by hitting Ctrl-C\n\n";
		sleep($safetime);
	}
		
	
	foreach my $pip (@piplst){
		undef $pipes{$pip};
	}
		
	pod2usage(-msg     => "\n\nNo directories to process\n\n",
		  -exitval => 1,
	          -verbose => 1,
		  -output  => \*STDERR) unless keys %pipes;
	
	
        if (defined $fragment_length &! $fragment_length){
                
                foreach my $type qw(aat_na aat_aa){
		        my $basedir = "$pipdir/$type";
			
		        print STDERR "Type: $type\nBasedir: $basedir\n" if $debug == 1;

		        next unless -d $basedir;

		        print STDERR "\tExists\n" if $debug > 1;

		        if (opendir(my $resdir, $basedir)){
			        
                                foreach my $sdir (readdir($resdir)){

				        print STDERR "Result Dir: $sdir\n" if $debug > 1;

				        next unless $sdir =~ /(\d+)_(\S+)/ && exists $pipes{$1};
				        my ($pip, $searchdb) = ($1, $2);
				        my $tfile = $searchdb;
					
					print STDERR "Result Dir: $basedir/$sdir\n" if $debug == 1;
					
                                        while (my ($outtype, $listname) = each %{$list_names{$type}}){
					
						print STDERR "Output type: $outtype\n" if $debug == 1;;
						
                                                my $list  = "$basedir/$sdir/$listname";
				                my $tfile = $outtype eq 'btab' ? "$searchdb.$tool_tr{$type}.btab" : "$searchdb.$tool_tr{$type}";

				                print STDERR "\tExists\nTarget File: $tfile\nList file: $list" if $debug > 1;

				                if (-s "$list" && -r "$list"){
					                print STDERR "\tExists\n" if $debug > 1;

					                my $cmd = "$mover -D $db -l 0 -L $list -t $tfile $move_opt{$outtype}";
        					        $cmd .= ' -r' if $remove;
	        				        print "CMD: '$cmd'\n";
        
	        				        system($cmd) && warn "Pipeline ID: '$pip' Result type: '$type'\t- Some error has occoured:\t\"$!\"\n\n";
		        		        } else {
			        		        print "Impossible to find the list '$list'\n\n";
				                        #warn "Impossible to find the list '$list'\n\n";
				                }
                                        }
			        }
			        closedir($resdir);
		        } else {
			        warn "Impossible to open the directory '$basedir' - skipping\n";
		        }
	        }
        } else {
                foreach my $type qw(btab gap2 nap){
		        my $basedir = "$pipdir/output_repository/adjust_${type}_coordinates";

		        print STDERR "Basedir: $basedir\n" if $debug == 1;

		        next unless -d $basedir;

		        print STDERR "\tExists\n" if $debug > 1;

		        if (opendir(my $resdir, $basedir)){
			        my $listfile = "adjust_${type}_coordinates.list";

			        print STDERR "List File Name: $listfile\n" if $debug == 1;

			        foreach my $sdir (readdir($resdir)){

				        print STDERR "Result Dir: $sdir\n" if $debug == 1;

				        next unless -d "$basedir/$sdir" && $sdir =~ /(\d+)_(\S+)/ && exists $pipes{$1};
				        my ($pip, $searchdb) = ($1, $2);
				        my $list  = "$basedir/$sdir/$listfile";
				        my $tfile = $searchdb;

				        print STDERR "\tExists\nTarget File: $tfile\nList file: $list" if $debug > 1;

				        if (-s "$list" && -r "$list"){
					        print STDERR "\tExists\n" if $debug > 1;

					        if ($type eq 'btab'){ # in the case of btab files, we need to open the list to figure out if they result from gap2 or nap searches
						        open(my $lst, $list) || die "\n\nProblems opening the file: '$list' for reading\n\n";
						        chomp(my $firstfile = <$lst>);
						        close($lst);
						        my ($stype) = ($firstfile =~ /(gap2|nap)/);
						        $tfile .= ".$stype.btab";
					        } else {
						        $tfile .= ".adj.$type";
					        }

					        my $cmd = defined $fragment_length ? "$mover -D $db -l $fragment_length -L $list -t $tfile $move_opt{$type}" : "$mover -D $db -L $list -t $tfile $move_opt{$type}";
                          		        $cmd .= ' -r' if $remove;
					        print "CMD: '$cmd'\n";

					        system($cmd) && warn "Pipeline ID: '$pip' Result type: '$type'\t- Some error has occoured:\t\"$!\"\n\n";
				        } else {
					        warn "Impossible to find the list '$list'\n\n";
				        }
			        }
			        closedir($resdir);

		        } else {
			        warn "Impossible to open the directory '$basedir' - skipping\n";
		        }
	        }
        }
}
