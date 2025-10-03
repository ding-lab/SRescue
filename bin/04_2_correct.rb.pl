#!/usr/bin/env perl
$workdir=$ARGV[0];
open(supread,"$workdir/supread.bam.tsv");
while(<supread>)
{
	chomp;
	@token=split(/\t/,$_);
	if(grep(/SA/,$_))
	{
		$id=$token[0];
		$chr1=$token[2];
		$start1=$token[3];

		@token1=split(/SA:Z:/,$_);
		@token2=split(/\;/,$token1[1]);
		@token3=split(/\,/,$token2[0]);
		$chr2=$token3[0];
		$start2=$token3[1];
		$strand2=$token3[2];
		# @token4=split(/^(\d+)/,$token3[3]);
		$frag2=$token3[3];
		$frag2 =~ /^(\d+)/;
		$frag2=$1;
		# print $token3[3]."\t".$1."\n";
		
		if($strand2 eq "+")
		{
			$end2=$start2;
		}else{
			$end2=$start2+$frag2;
		}
		if($chr1 eq $chr2)
		{
			if($start1<$start2)
			{
				# print $id."\t".$start1."\t".$start2."\t".$end2."\t".$frag2."\n";
				$nontra{$id}=$start2."\t".$frag2."\t".$strand2;
				${$id}{$chr2}=0;
				${$id}{$chr1}=0;
			}
		}else{
			${$id}{$chr2}=0;
			${$id}{$chr1}=0;
		}
		# else{
		# 	$tra{$id}=0;
		# 	${$id}{$chr2."\t".$end2}=0;
		# }
	}
	

}
close(file);



open(file,"$workdir/finalfilter_cutesv.vcf");
open(res,">$workdir/finalfilter_cutesv2.vcf");
while(<file>)
{
	chomp;
	# print $_."\n";
	@ss=split(/\t/,$_);
	if(grep(/^chr/,$_))
	{
		@token=split(/RNAMES=/,$_);
		@token1=split(/\;/,$token[1]);
		# print $token[1]."\n";
		@token1=split(/\,/,$token1[0]);
		if(grep(/SVTYPE=DEL/,$_) or grep(/SVTYPE=DUP/,$_))
		{
			$readid=$token1[0];
			if(exists($nontra{$readid}))
			{
				@getend=split(/END=/,$_);
				@getend1=split(/\;/,$getend[1]);
				$endpos=$getend1[0];

				@keys=keys(%{$readid});
				$num=scalar @keys;
				# print $readid."\t".$num."\n";
				if($num >= 2)
				{
					next;
				}else{
					$combo=$nontra{$readid};
					@scombo=split(/\t/,$combo);
					$start2=$scombo[0];
					$frag2=$scombo[1];
					$strand2=$scombo[2];
					if(grep(/SVTYPE=DEL/,$_))
					{
						if($strand2 eq "+")
						{
							$end=$start2;
						}else{
							$end=$start2;
						}
						$svlen=$ss[1]-$end;
					}else{
						if($strand2 eq "+")
						{
							$end=$start2+$frag2;
							# print $ss[2]."\t".$frag2."\n";
						}else{
							$end=$start2+$frag2;
							
						}
						# print $ss[1]."\t".$ss[2]."\t".$end."\t".$strand2."\n";
						$svlen=abs($ss[1]-$end);
					}
					
					@news0=split(/SVLEN=/,$_);
					@news1=split(/\;END=/,$news0[1]);
					@news2=split(/\;RE/,$news1[1]);
					print $ss[0]."\t".$ss[1]."\t".$end."\t".$endpos."\n";
					if($end> ($endpos-1000) and ($end<($endpos+1000)))
					{

						print res $news0[0]."SVLEN=".$svlen.";END=".$end.";RE".$news2[1]."\n";
					}

				}
			}else{
				print res $_."\n";
			}
		}else{
			print res $_."\n";
		}


	}else{
		print res $_."\n";
	}
}
close(file);
close(res);

open(file,"$workdir/supread.depth.tsv");
while(<file>)
{
	chomp;
	@token=split(/\t/,$_);
	$depth{$token[0]."\t".$token[1]}=$token[2];
	# print $_."\n";
}
close(file);

open(file,"$workdir/finalfilter_cutesv2.vcf");
open(res,">$workdir/finalfilter_cutesv3.vcf");
while(<file>)
{
	chomp;
	@token=split(/\t/,$_);
	if(grep(/^chr/,$_))
	{
		$ed=$token[1]+1;
		# print $token[0]."\t".$ed."\n";
		if(exists($depth{$token[0]."\t".$ed}))
		{
			$dp=$depth{$token[0]."\t".$ed};
			@token1=split(/\:/,$token[9]);
			$dv=$token1[2];
			$dr=$dp-$dv;
			$out9=$token1[0].":".$dr.":".$dv.":".$token1[3].":".$token1[4];
			$vaf=$dv/$dp;
			if(grep(/AF=/,$token[7]))
			{
				@token2=split(/AF=/,$token[7]);
				$out7=$token2[0]."AF=".$vaf;
			}else{
				$out7=$token[7].";".$vaf;
			}
			print res $token[0]."\t".$token[1]."\t".$token[2]."\t".$token[3]."\t".$token[4]."\t".$token[5]."\t".$token[6]."\t".$out7."\t".$token[8]."\t".$out9."\n";
		}
	}else{
		print res $_."\n";
	}
}
close(file);
close(res);
# perl /storage1/fs1/dinglab/Active/Projects/yuweiz/projects/HTAN/ccRCC/Longread/t2t/01_compare2sr/00_start_allSample/script/04_2_correct.rb.pl /storage1/fs1/dinglab/Active/Projects/yuweiz/projects/HTAN/ccRCC/Longread/t2t/01_compare2sr/00_start_allSample/HT941/
