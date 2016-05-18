$head=$ARGV[0];

open(IN ,"$head.bim") || die;
open(OUT,">$head.list") || die;
while (<IN>) {
	chomp $_; @line=split;
	if(defined($rs{$line[1]})) {
		printf OUT "$line[1]\n";
	} else {
		$rs{$line[1]}=1;
	}
}
close IN;
close OUT;