#!/usr/bin/env perl 

# VEP Slicer: Extract variants from Variant Effect Predictor (VEP) annotated VCFs by keyword
#
# Author: R. Jay Mashl <rjmashl@gmail.com>
# 2016-06-16: v0.1: initial version
#
# Inputs: VCF file and one or more VEP keywords
# Output: Variant calls (chr,pos,ref,alt) with VEP values for the keywords

use strict;
use warnings;

my $vcf;          # VCF file
my @label = ();     # list of VEP keywords
my @fixed_col = (); # 0-based column index of VEP keywords 

sub get_value {
    my $info_str   = shift;
    my $field_name = shift;
    
    my @info_arr = split /\;/,$info_str;
    for(my $i = 0; $i < scalar @info_arr; $i++) {
	if( index($info_arr[$i], $field_name) == 0) { # match start of string
	    my @b = split /=/, $info_arr[$i];
	    return $b[1];
	}
    }
    return "NULL";  # no match
}

sub syntax {
    my $sn = `basename $0`;
    chomp( $sn );
    my $usage = qq{     Usage:

     $sn   vcffile

        - display list of VEP keywords, if present, from vcffile

     $sn   vcffile  VEP_keyword  [VEP_keyword] ...

        - display variant calls (chr,pos,ref,alt) from vcffile and VEP values for the keywords
};

    die $usage;

}

sub get_vep_keywords {
    my $vcf = shift;
    print "\n";
    open(IN , "< $vcf")  ||  die "\nError: Cannot find file $vcf\n\n\n";
    while(<IN>) {
	chomp;
        if(/^##INFO=<ID=CSQ/) {
            my @info_id = split /Format:\ /, $_;
            $info_id[1] = substr( $info_id[1], 0, -2 );  # remove end xml chars
            my @tmp     = split /\|/, $info_id[1];
	    print "VEP Keywords present:\n";
	    print join(", ", @tmp);
	}
    }
    close(IN);
    print "\n\n";
}



if( scalar @ARGV == 0) {
    syntax;
}
if( scalar @ARGV == 1) {
    get_vep_keywords $ARGV[0];
    syntax;
}


($vcf, @label) = @ARGV;
$fixed_col[0] = -1;  # flag denoting column index is unassigned


open(IN , "< $vcf")  ||  die "\nError: Cannot find file $vcf\n\n\n";
while(<IN>) {
    chomp;
    if( $fixed_col[0] < 0 ) {   # check if have not seen VEP info
	if(/^##INFO=<ID=CSQ/) {
	    # get field map, e.g.
            ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|ASN_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE">
	    
	    my @info_id = split /Format:\ /, $_;
	    $info_id[1] = substr( $info_id[1], 0, -2 );  # remove end xml chars
	    my @tmp     = split /\|/, $info_id[1];
	    
	    # searching/setting column in input order
	    my $header_str = "#CHROM\tPOS\tREF\tALT";
	    for( my $j = 0; $j < scalar @label ;  $j++) {
		for( my $i = 0; $i < scalar @tmp ;  $i++) {
		    if( $tmp[$i] eq $label[$j] ) {
			$fixed_col[$j]  = $i;
			$header_str    .= "\t".$tmp[$i];
		    }
		}
	    }
	    print $header_str, "\n";
	}
	else {
	    next;
	}
    }

    # skip subsequent comments, if any
    next if( /^#/ );   

    # process calls
    my @field = split /\t/, $_;
    my $value = get_value( $field[7], "CSQ" );
    
    if( $value ne "NULL" ) {  # CSQ was present
	my @anno_field = split /\|/, $value, -1;
	
	my $print_str = join( "\t", @field[0,1,3,4] );
	for(my $j = 0 ; $j < scalar @fixed_col ; $j++) {
	    my $answer = $anno_field[ $fixed_col[$j] ];
	    if( length($answer) == 0 ) {
		$answer = "NULL";
	    }
	    $print_str .= "\t".$answer;
	}
	print $print_str, "\n";
    }
}
close( IN );
