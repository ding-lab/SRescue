open(normal,"/storage1/fs1/dinglab/Active/Projects/yuweiz/projects/HTAN/ccRCC/Longread/t2t/01_compare2sr/00_start_allSample/HT941/LRcaller_normal.call_noins.vcf");
while(<normal>)
{
	chomp;
	@token=split(/\t/,$_);
	if(grep(/^chr/,$_))
	{
		@token1=split(/\:/,$token[-1]);
		@token2=split(/\,/,$token1[-2]);
		print $token2[1]."\n";
		if($token2[1] == 0)
		{
			$normal{$token[2]}=0;
		}
		
	}
}
close(normal);

open(file,"/storage1/fs1/dinglab/Active/Projects/yuweiz/projects/HTAN/ccRCC/Longread/t2t/01_compare2sr/00_start_allSample/HT941/LRcaller_tumor.call_noins.vcf");
open(res,">/storage1/fs1/dinglab/Active/Projects/yuweiz/projects/HTAN/ccRCC/Longread/t2t/01_compare2sr/00_start_allSample/HT941/finalfilter.LRCaller.vcf");
while(<file>)
{
	chomp;
	@token=split(/\t/,$_);
	if(grep(/^chr/,$_))
	{
		@token1=split(/\:/,$token[-1]);
		@token2=split(/\,/,$token1[-2]);
		$alt=$token2[1];
		$total=$token2[2];
		print $token[2]."\t".$alt."\t".$total."\n";
		if($total>0)
		{
			$vaf=$alt/$total;
			if(exists($normal{$token[2]}) and ($alt>0))
			{
				print res $token[0]."\t".$token[1]."\t".$token[2]."\t".$token[3]."\t".$token[4]."\t".$token[5]."\t".$token[6]."\t".$token[7]."\t".$token[8]."\t".$token[15]."\n";
			}
		}

	}else{
		print res $_."\n";
	}
}
close(res);
close(file);
