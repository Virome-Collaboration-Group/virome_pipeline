#!/usr/bin/perl -w

=head1 NAME

mga2seq_pep.pl - convert metagene raw output to fasta nucliotide and peptide file

=head1 SYNOPSIS

USAGE: mga2seq_pep.pl
            --input=/path/to/fasta
            --mga=/path/to/mga/output
    	    --prefix
    	    --outdir
          [ --log=/path/to/logfile
            --debug=N
          ]

=head1 OPTIONS

B<--input, -i>
    The full path to fasta sequence file.

B<--input, mga>
    Metagene output file

B<--debug,-d>
    Debug level.  Use a large number to turn on verbose debugging.

B<--log,-l>
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

    This script is used to convert metagene output to nucleotide and peptide seq file.

=head1  INPUT

    The input to this is defined using the --input/-i or mga/-m.  This should point
    to the fasta file containing sequence(s), and metagene output file

=head1 AUTHOR

    Written by Daniel Nasko,
    Center for Bioinformatics and Computational Biology, University of Delaware.

=head1 COPYRIGHT

    Copyright 2017 Daniel Nasko.
    License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.
    This is free software: you are free to change and redistribute it.
    There is NO WARRANTY, to the extent permitted by law.

    Please acknowledge author and affiliation in published work arising from this script's
    usage <http://bioinformatics.udel.edu/Core/Acknowledge>.

=cut

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

BEGIN {
	use Ergatis::Logger;
}

my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
						  'mga|m=s',
                          'prefix|p=s',
						  'outdir|o=s',
                          'log|l=s',
                          'debug|d=s',
                          'help|h') || pod2usage();

#### display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

##############################################################################
#### make sure everything passed was peachy
&check_parameters(\%options);

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

my %Fasta;
my ($h,$s,$l) = ("","",0);
my %Model = (
    'a' => 'archaea',
    'p' => 'phage',
    'b' => 'bacteria',
    's' => 'self');
my %Codon;
my $nseqs=0;
LoadCodonTable();

open(IN, "<", $options{input}) || die "\n\n Cannot open the input file: $options{input}\n\n";
while(<IN>) {
    chomp;
    if ($_ =~ m/^>/) {
	unless ($l == 0) {
	    $s = uc($s);
	    $Fasta{$h} = $s;
	    $s = "";
	}
	$h = $_;
	$h =~ s/^>//;
	$h =~ s/ .*//;
    }
    else {
	$s = $s . $_;
    }
    $l++;
}
close(IN);
$s = uc($s);
$Fasta{$h} = $s;

#### get metagene output file based on prefix
my $file = `egrep "$options{prefix}\\." $options{mga}`;   # modified 4/4/11 SWP

unless (length($file)){
    $logger->logdie("No metagene output file for $options{prefix}");
    exit(1);
}

# NODE_1_length_56803_cov_57.912_g0_i0_420_4775_1 size=1452 gc=0.302607 start=420 stop=4775 strand=+ frame=0 model=self score=844.516 type=complete caller=MetaGENE

open(PEP, ">" , "$options{outdir}/$prefix.pep") || $logger->logdie("Cannot write to: $options{outdir}/$prefix.pep");
open(NUC, ">", "$options{outdir}/$prefix.nuc") || $logger->logdie("Cannot write to: $options{outdir}/$prefix.nuc");

open(IN, "<", $file) || $logger->logdie("Cannot open the mga file: $file");
while(<IN>) {
    chomp;
    if ($_ =~ m/^gene/) {
    	my @a = split(/\t/, $_);
    	my $gid = get_gid($a[0]);

        if (exists $Fasta{$h}) {
    	    my $nt_orf_seq = get_nt_orf($Fasta{$h}, $a[1], $a[2], $a[3], $a[4]);
    	    if ($a[3] eq "-") { my $tmp=$a[1]; $a[1]=$a[2];$a[2]=$tmp; }
    	    print NUC ">" . $h . "_" . $a[1] . "_" . $a[2] . "_" . $gid . " size=" . length($nt_orf_seq) . " gc=" . gc($nt_orf_seq) . " start=$a[1]" . " stop=$a[2]" . " strand=$a[3]" . " frame=$a[4]" . " model=" . $Model{$a[7]} . " score=$a[6]" . " type=" . get_type($a[5]) . " caller=MetaGENE" . "\n";
    	    print NUC $nt_orf_seq . "\n";

    	    print PEP ">" . $h . "_" . $a[1] . "_" . $a[2] . "_" . $gid . " size=" . orf_len($nt_orf_seq) . " gc=" . gc($nt_orf_seq) . " start=$a[1]" . " stop=$a[2]" . " strand=$a[3]" . " frame=$a[4]" . " model=" . $Model{$a[7]} . " score=$a[6]" . " type=" . get_type($a[5]) . " caller=MetaGENE" . "\n";
    	    my $translated_seq = translate($nt_orf_seq);
    	    if ($a[5] =~ m/^1/) {$translated_seq =~ s/^./M/;} ## This seems weird, but isnt. MetaGene allows for alternative start codons when predicting ORFs, funny thing about alternative start codons is, even if they are supposed to encode some other peptide the tRNA will always put a Met there.
    	    print PEP $translated_seq . "\n";
    	    $nseqs++;
	    } else {
            $logger->logdie("Error! Cannot find the sequence: $h");
        }
    } elsif ($_ =~ m/^#/) {
        unless ($_ =~ m/^# gc = / || $_ =~ m/^# self: /) {
            $h = $_;
            $h =~ s/^# //;
            $h =~ s/ .*//;
        }
    }
}
close(IN);
close(PEP);
close(NUC);

if ($nseqs == 0) { die "\n Error: There were no ORFs predicted from your input files\n FASTA file: $fasta\n MGA file: $mga\n\n"; }

##############################################################################
sub translate
{
    my $str = $_[0];
    my $peptide = "";
    for (my $i=0; $i<length($str)-2; $i+=3) {
    	my $mer = substr $str, $i, 3;
    	if (exists $Codon{$mer}) {
            $peptide = $peptide . $Codon{$mer};
        } else {
            $peptide = $peptide . 'X';
        }
    }

    return $peptide;
}

##############################################################################
sub LoadCodonTable
{   # Load the AA codon table from __END__ of program.
    my @data = <DATA>;
    foreach my $line (@data) {
        chomp($line);
	    my @codons = split(/ /,$line);
	    my $AA = shift(@codons);
	    foreach my $nnn (@codons) {
            $nnn =~ s/U/T/g;
		    $Codon{$nnn} = $AA;
	    }
    }
}

##############################################################################
sub revcomp
{
    my $str = $_[0];
    $str = scalar reverse $str;
    $str =~ tr/ATGC/TACG/;
    return $str;
}

##############################################################################
sub get_type {
    my $str = $_[0];

    if ($str eq "11") {
        return "complete";
    }
    elsif ($str eq "10") {
        return "lack_stop";
    }
    elsif ($str eq "01") {
        return "lack_start";
    }
    elsif ($str eq "00") {
        return "incomplete";
    }
}

##############################################################################
sub gc {
    my $str = $_[0];
    my $gcs = $str =~ tr/GCgc/GCGC/;
    my $gc_content = $gcs / length($str);
    return $gc_content;
}

##############################################################################
sub orf_len {
    my $str = $_[0];
    my $len = length($str)/3;
    return $len;
}

##############################################################################
sub get_nt_orf {
    my $seq = $_[0];
    my $start = $_[1];
    my $stop = $_[2];
    my $sense = $_[3];
    my $frame = $_[4];

    my $length = $stop - $start + 1;
    $start--;

    my $orf = substr $seq, $start, $length;

    if ($sense eq "-") {
        $orf = revcomp($orf);
    }
    $orf = substr $orf, $frame;

    return $orf;
}

##############################################################################
sub get_gid {
    my $str = $_[0];
    $str =~ s/gene_//;
    return $str;
}

##############################################################################
sub check_parameters {
    #### make sure sample_file and output_dir were passed
    my @required = qw(input mga outdir);

	foreach my $key (@required) {
		unless ($options{$key}) {
			pod2usage({-exitval => 2,  -message => "ERROR: Input $key not defined", -verbose => 1, -output => \*STDERR});
		    $logger->logdie("No input defined, plesae read perldoc $0\n\n");
        }
	}
}

__END__
A GCT GCC GCA GCG
R CGT CGC CGA CGG AGA AGG
N AAT AAC
D GAT GAC
C TGT TGC
Q CAA CAG
E GAA GAG
G GGT GGC GGA GGG
H CAT CAC
I ATT ATC ATA
L TTA TTG CTT CTC CTA CTG
K AAA AAG
M ATG
F TTT TTC
P CCT CCC CCA CCG
S TCT TCC TCA TCG AGT AGC
T ACT ACC ACA ACG
W TGG
Y TAT TAC
V GTT GTC GTA GTG
* TAA TAG TGA
