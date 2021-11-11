#!/usr/bin/perl
use Switch;
use Time::HiRes qw( time );
my $PDBID=$ARGV[0];
my $deschain=$ARGV[1];
my $tgtchain=$ARGV[2];
my $LIBPDBID=$ARGV[3];
#####PARAMETERS FROM LAUNCH INTERFACE########
my $dist_cutoff=4.5; #"4.5";
my $ialign_type='is'; #"is";

my $mini='15'; #"15";
#############################################



my $fh_des_resids="$PDBID.$deschain"; 
my $fh_tgt_resids="$PDBID.$tgtchain"; 
my $IALIGN="ialign/bin/ialign.pl"; 
my $ial_lib_rescol=7; #column of the amino acid code of the library complex in ialign
my $ial_tgt_rescol=3; #column of the amino acid code of the library complex in ialign

%des_interface=();
%tgt_interface=();


my $des_res_file="$PDBID"."_"."$deschain".'.res';
my $tgt_res_file="$PDBID"."_"."$tgtchain".'.res';
my $start = time();
`echo "$LIBPDBID">>pdb.lst`;
my $ial_msa="$PDBID"."$deschain".'.ial_msa';
my $ial_tarmsa="$PDBID"."$tgtchain".'.ial_msa';
open(IALMSA, ">>$ial_msa");
open(IALTARMSA, ">>$ial_tarmsa");



`$IALIGN  -a 2 -minp 4 -dc $dist_cutoff -e $ialign_type -mini $mini $PDBID.pdb $LIBPDBID.pdb >ial_temp.txt`;
#exit early if a valid interface is not found
$no_PPI=`grep -c 'found 0 valid PPI' 'ial_temp.txt'`;

if ($no_PPI ==1) {
	`rm ${LIBPDBID}*`;
	exit;
}
print IALMSA "$LIBPDBID\t";
print IALTARMSA "$LIBPDBID\t";

open(IALTEMP, "<ial_temp.txt");
while($line = <IALTEMP>) {
	if( ($line =~ /IS-score/) || ($line =~ /TM-score/) ) {
		chomp(my $iscore=(split /[\s+,]/, $line)[2]);
		print IALMSA "$iscore\t";
		print IALTARMSA "$iscore\t";
	}

	if ($line=~m/\s+\w+\s+\w{1}\s+\d+\s+\w{3}\s+\w\s+\d+\s+\w{3}\s+\d+.\d+\s+\d+/ ){
		@align=split /\s+/, $line;		
		my $aa=changeAA($line,$ial_lib_rescol);
		if ($line=~m/\s+\w+\s+$tgtchain\s+\d+\s+\w{3}\s+\w\s+\d+\s+\w{3}\s+\d+.\d+\s+\d+/ ){
			$lib_tgtchain=@align[5];
			$tgt_interface{@align[$ial_tgt_rescol]}=$aa;
		}
		if ($line=~m/\s+\w+\s+$deschain\s+\d+\s+\w{3}\s+\w\s+\d+\s+\w{3}\s+\d+.\d+\s+\d+/ ){
			$lib_deschain=@align[5];
			$des_interface{@align[$ial_tgt_rescol]}=$aa;
		}
	}
}

close IALTEMP;
close IALMSA;
close IALTARMSA;	
print_int_al('target',%tgt_interface);
print_int_al('design',%des_interface);


`rm ${LIBPDBID}.pdb`;
`rm *con.lst`;
`rm *int.pdb`;
`rm ${LIBPDBID}.parsed`;

sub changeAA {
	$line=@_[0];
	$aa_pos=@_[1];
#	print "@_[1] @_[0]";
	@word_array=split /\s+/, $line;
#	print "achange  $aa_pos @word_array[$aa_pos]\n";
	switch (@word_array[$aa_pos]){
		case "ALA"  {return 'A'}
		case "ARG"  {return 'R'}
		case "ASP"  {return 'D'}
		case "ASN"  {return 'N'}
		case "CYS"  {return 'C'}
		case "GLN"  {return 'Q'}				
		case "GLY"  {return 'G'}
		case "GLU"  {return 'E'}
		case "HIS"  {return 'H'}
		case "ILE"  {return 'I'}
		case "LEU"  {return 'L'}
		case "LYS"  {return 'K'}
		case "MET"  {return 'M'}
		case "PHE"  {return 'F'}
		case "PRO"  {return 'P'}
		case "SER"  {return 'S'}
		case "THR"  {return 'T'}
		case "TRP"  {return 'W'}
		case "TYR"  {return 'Y'}
		case "VAL"  {return 'V'}	
	}
}






sub print_int_al() {
	my $chain=@_[0];

	my  %interface=@_[1];
	if ($chain eq 'design') {
		my $res_file=$des_res_file;
		my $msa_file=$ial_msa;
		open RES, "<$res_file";
		chomp(my @resids = <RES>);
		open FHMSA, ">>$msa_file";
		foreach ( @resids) {	
			if (exists $des_interface{$_} ) { 
				print FHMSA $des_interface{$_};
			} 
			else {	
				print FHMSA "-";
			}
		}
	print FHMSA "\n";
	}
	if ($chain eq 'target') {
		my $res_file=$tgt_res_file;
		my $msa_file=$ial_tarmsa;
		open RES, "<$res_file";
		chomp(my @resids = <RES>);
		open FHMSA, ">>$msa_file";
		foreach ( @resids) {	
			if (exists $tgt_interface{$_} ) { 
				print FHMSA $tgt_interface{$_};
			} 
			else {	
				print FHMSA "-";
			}
		}
	print FHMSA "\n";
	}
	close FHMSA;
	close RES;
}







