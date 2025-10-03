#!/usr/bin/env perl
$input=$ARGV[0];
open(file,"$input");
$output=$input;
$output=~s/\.vcf/4fc.vcf/g;
open(res,">$output");
while(<file>)
{
	chomp;
	@token=split(/\t/,$_);
	if(!grep(/^chr/,$_))
	{
		print res $_."\n";
	}else{
		if(grep(/SVTYPE=DUP/,$_))
		{
			print res $_."\n";
			$line=$_;
			$line=~s/SVTYPE=DUP/SVTYPE=INS/g;
			$line=~s/<DUP>/<INS>/g;
			print res $line."\n";
		}else{
			print res $_."\n";
		}
	}
}
close(res);
close(file);
