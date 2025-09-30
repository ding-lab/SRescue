$path=$ARGV[0];
$tumor=$path."/cutesv_tumor.call.vcf";
$normal=$path."/cutesv_normal.call.vcf";
$outfile=$path."/finalfilter_cutesv.vcf";
$supread=$path."/supread.tsv";
$loci=$path."/supread.loci.bed";

open(normal,"$normal");
while(<normal>)
{
	chomp;
	if(grep(/^chr/,$_))
	{
		@token=split(/\t/,$_);
		@token1=split(/\:/,$token[9]);
		$dv=$token1[2];
		if($token[4] eq "<DUP>")
		{
			$hash{$token[0]."\t".$token[1]}=$dv;
		}elsif($token[4] eq "<INS>")
		{
			if(exists($hash{$token[0]."\t".$token[1]}))
			{
				$dupdv=$hash{$token[0]."\t".$token[1]};
				if($dv>$dupdv)
				{
					$hash{$token[0]."\t".$token[1]}=$dv;
				}

			}else{
				$hash{$token[0]."\t".$token[1]}=$dv;
			}
			# $hash{$token[0]."\t".$token[1]}=0;
		}else{
			# print $_."\t".$dv."\n";
			$hash{$token[0]."\t".$token[1]}=$dv;
		}
	}
}
close(normal);


open(res,">$outfile");
open(tumor,"$tumor");
while(<tumor>)
{
	chomp;
	if(grep(/^chr/,$_))
	{
		@token=split(/\t/,$_);
		@token1=split(/\:/,$token[9]);
		$dv=$token1[2];
		$normaldv=$hash{$token[0]."\t".$token[1]};
		if($normaldv == 0)
		{
			if($token[4] eq "<DUP>")
			{
				$dup{$token[0]."\t".$token[1]}=$_;
				# print $token[0]."\t".$token[1]."\n";
			}
			if($token[4] eq "<INS>")
			{
				$ins{$token[0]."\t".$token[1]}=$_;
				# if(exists($tumor{$token[0]."\t".$token[1]}))
				# {
				# print $token[0]."\t".$token[1]."\t".$dv."\n";

				# 	@token2=split(/\:/,$tumor{$token[0]."\t".$token[1]});
				# 	$dupdv=$token2[-1];
				# 	if($dv>$dupdv)
				# 	{
				# 		$out=$_;
				# 		$out=~s/INS/DUP/g;
				# 		$tumor{$token[0]."\t".$token[1]}=$out;
				# 	}

				# }else{
				# 	$tumor{$token[0]."\t".$token[1]}=$_;
				# }
			}
			if($token[4] ne "<DUP>" and $token[4] ne "<INS>"){
				# print $_."\n";
				$other{$token[0]."\t".$token[1]}=$_;
			}
		}
	}else{
		print res $_."\n";
	}
}
close(tumor);
while(my ($k,$v)=each %hash)
{
	# @tt=split(/\:/,$other{$k});
	$normaldv=$v;
	# print $normaldv."\n";
	if(exists($other{$k}))
	{
		@ss=split(/\t/,$other{$k});
		@token2=split(/\:/,$ss[9]);
		$otherdv=$token2[2];

		# print $other{$k}."\n";
		if($otherdv > 0)
		{
			# print $k."\t".$otherdv."\n";
			print res $other{$k}."\n";
		}
	}
	if(exists($ins{$k}))
	{
		@ss=split(/\t/,$ins{$k});
		@token2=split(/\:/,$ss[9]);
		$insdv=$token2[2];
		# print $ins{$k}."\n";
		$output=$ins{$k};
		if(exists($dup{$k}))
		{
			@ss=split(/\t/,$dup{$k});
			@token1=split(/\:/,$ss[9]);
			$dupdv=$token1[2];
			# print $insdv."\t".$dupdv."\n";
			if($insdv>$dupdv and $insdv>0)
			{
				# print $output."\n";
				$output=~s/INS/DUP/g;
				print res $output."\n";
			}else{
				# print $dup{$k}."\n";				
				if($dupdv > 0)
				{
					print res $dup{$k}."\n";
				}
			}
		}else{
			if($dupdv > 0)
			{
			print res $ins{$k}."\n";
		}
		}
	}
	
}
close(res);
# print $outfile."\n";
open(file,"$outfile");
open(res,">$supread");
open(res1,">$loci");
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
		print res $token1[0]."\n";
		$send=$ss[1]+1;
		print res1 $ss[0]."\t".$ss[1]."\t".$send."\n";

	}
}
close(file);
close(res);

