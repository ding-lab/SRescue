$INPUT=$ARGV[0];
$OUTPUT=$ARGV[1];
$OUTPUT2=$ARGV[2];

open(file,"$INPUT");
open(res,">$OUTPUT");
open(res1,">$OUTPUT2");

%hash=();
$n=0;

while(<file>)
{
	chomp;
	@token=split(/\t/,$_);
	if(grep(/^chr/,$_))
	{
		$chr1=$token[0];
		$loc1=$token[1];
		@token1=split(/END=/,$_);
		@token2=split(/\;/,$token1[1]);
		$loc2=$token2[0];
		@token3=split(/CHR2=/,$_);
		@token4=split(/\;/,$token3[1]);
		$chr2=$token4[0];
		@token5=split(/SVTYPE=/,$_);
		@token6=split(/\;/,$token5[1]);
		$type=$token6[0];
		@token7=split(/SVLEN=/,$_);
		@token8=split(/\;/,$token7[1]);
		$len=abs($token8[0]);
		$id="SURVIVOR"."_".$n;
		$n++;
		$flag="T";
		if($token[9] eq "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN")
		{
			$sr_flag="FALSE";
		}else{
			@sv1=split(/\:/,$token[9]);
			$sr_flag=$sv1[7];
			$sr_type=$sv1[6];
		}
		if($token[10] eq "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN")
		{
			$lr_flag="FALSE";
		}else{
			@sv2=split(/\:/,$token[10]);
			$lr_flag=$sv2[7];
			$lr_type=$sv2[6];
		}

		if(grep(/\>/,$token[4]))
		{
			$ts=$token[4];
			$ts=~s/\>//g;
			$ts=~s/\<//g;
			$type=$ts;
			# print $type."\n";
		}
		# $loc1_1=$loc1-1;
		# $loc1_2=$loc1+1;
		# $loc2_1=$loc2-1;
		# $loc2_2=$loc2+1;
		if($type ne "TRA")
		{
			if($type eq "INV")
			{
				$svlen=$loc2-$loc1;
				$len="SVLEN=".$svlen;
				# print $len."\n";
				$line=$_;
				$line=~s/SVLEN=0/$len/g;
				# print $line."\n";
				@token=split(/\t/,$line);
				@t1=split(/SVTYPE=/,$token[7]);
				@t2=split(/\;SVME/,$t1[1]);
				# print $t2[1]."\n";
				$resline=$t1[0]."SVTYPE=".$type.";SVME".$t2[1];
				print res1 $token[0]."\t".$token[1]."\t".$token[2]."\t".$token[3]."\t"."<".$type.">"."\t".$token[5]."\t".$token[6]."\t".$resline."\t".$token[8]."\t".$token[9]."\t".$token[10]."\n";
			}else{
				@t1=split(/SVTYPE=/,$token[7]);
				@t2=split(/\;SVME/,$t1[1]);
				$resline=$t1[0]."SVTYPE=".$type.";SVME".$t2[1];
				print res1 $token[0]."\t".$token[1]."\t".$token[2]."\t".$token[3]."\t"."<".$type.">"."\t".$token[5]."\t".$token[6]."\t".$resline."\t".$token[8]."\t".$token[9]."\t".$token[10]."\n";
			}
		}else{
			@t1=split(/SVTYPE=/,$token[7]);
			@t2=split(/\;SVME/,$t1[1]);
			$resline=$t1[0]."SVTYPE=".$type.";SVME".$t2[1];
			print res1 $token[0]."\t".$token[1]."\t".$token[2]."\t".$token[3]."\t".$token[4]."\t".$token[5]."\t".$token[6]."\t".$resline."\t".$token[8]."\t".$token[9]."\t".$token[10]."\n";
			# print res1 $token[0]."\t".$token[1]."\t".$token[2]."\t".$token[3]."\t"."<".$type.">"."\t".$token[5]."\t".$token[6]."\t".$resline."\t".$token[8]."\t".$token[9]."\t".$token[10]."\n";
		}

		if($type eq "INV")
		{
			$len=abs($loc2-$loc1);
		}

		# if(!exists($hash{$chr1."\t".$loc1}) and !exists($hash{$chr1."\t".$loc1_1}) and !exists($hash{$chr1."\t".$loc1_2}) and !exists($hash{$chr2."\t".$loc2}) and !exists($hash{$chr2."\t".$loc2_1}) and !exists($hash{$chr2."\t".$loc2_2}))
		if(!exists($hash{$chr1."\t".$loc1})  and !exists($hash{$chr2."\t".$loc2}))
		{
			@token_1=split(/STRANDS=/,$_);
			@token_2=split(/\;/,$token_1[1]);
			@token_3=split(//,$token_2[0]);
			$strand1=$token_3[0];
			$strand2=$token_3[1];
			$hash{$chr1."\t".$loc1}=$token[2];
			$hash{$chr2."\t".$loc2}=$token[2];
			# $hash{$chr1."\t".$loc1_1}=$token[2];
			# $hash{$chr1."\t".$loc1_2}=$token[2];
			# $hash{$chr2."\t".$loc2_1}=$token[2];
			# $hash{$chr2."\t".$loc2_1}=$token[2];
			$newhash{$token[2]}=$chr1."\t".$loc1."\t".$loc1."\t".$chr2."\t".$loc2."\t".$loc2."\t".$type."\t".$id."\t".$len."\t".$sr_flag."\t".$lr_flag."\t".$strand1."\t".$strand2;

		}else{
			if($type eq "INV")
			{
				if(exists($hash{$chr1."\t".$loc1}) and  grep(/SUPP=2/,$_))
				{
					$newk=$hash{$chr1."\t".$loc1};
					delete $newhash{$newk};
					$newhash{$token[2]}=$chr1."\t".$loc1."\t".$loc1."\t".$chr2."\t".$loc2."\t".$loc2."\t".$type."\t".$id."\t".$len."\t".$sr_flag."\t".$lr_flag."\t".$strand1."\t".$strand2;
				}
				# if(exists($hash{$chr1."\t".$loc1_1}) and  grep(/SUPP=2/,$_))
				# {
				# 	$newk=$hash{$chr1."\t".$loc1_1};
				# 	delete $newhash{$newk};
				# 	$newhash{$token[2]}=$chr1."\t".$loc1."\t".$loc1."\t".$chr2."\t".$loc2."\t".$loc2."\t".$type."\t".$id."\t".$len."\t".$sr_flag."\t".$lr_flag."\t".$strand1."\t".$strand2;
				# }
				# if(exists($hash{$chr1."\t".$loc1_2}) and  grep(/SUPP=2/,$_))
				# {
				# 	$newk=$hash{$chr1."\t".$loc1_2};
				# 	delete $newhash{$newk};
				# 	$newhash{$token[2]}=$chr1."\t".$loc1."\t".$loc1."\t".$chr2."\t".$loc2."\t".$loc2."\t".$type."\t".$id."\t".$len."\t".$sr_flag."\t".$lr_flag."\t".$strand1."\t".$strand2;
				# }
				if(exists($hash{$chr2."\t".$chr2}) and  grep(/SUPP=2/,$_))
				{
					$newk=$hash{$chr2."\t".$chr2};
					delete $newhash{$newk};
					$newhash{$token[2]}=$chr1."\t".$loc1."\t".$loc1."\t".$chr2."\t".$loc2."\t".$loc2."\t".$type."\t".$id."\t".$len."\t".$sr_flag."\t".$lr_flag."\t".$strand1."\t".$strand2;
				}
				# if(exists($hash{$chr2."\t".$loc2_1}) and grep(/SUPP=2/,$_))
				# {
				# 	$newk=$hash{$chr2."\t".$loc2_1};
				# 	delete $newhash{$newk};
				# 	$newhash{$token[2]}=$chr1."\t".$loc1."\t".$loc1."\t".$chr2."\t".$loc2."\t".$loc2."\t".$type."\t".$id."\t".$len."\t".$sr_flag."\t".$lr_flag."\t".$strand1."\t".$strand2;
				# }
				# if(exists($hash{$chr2."\t".$loc2_2}) and grep(/SUPP=2/,$_))
				# {
				# 	$newk=$hash{$chr2."\t".$loc2_2};
				# 	delete $newhash{$newk};
				# 	$newhash{$token[2]}=$chr1."\t".$loc1."\t".$loc1."\t".$chr2."\t".$loc2."\t".$loc2."\t".$type."\t".$id."\t".$len."\t".$sr_flag."\t".$lr_flag."\t".$strand1."\t".$strand2;
				# }
			}
		}


		# }
	}
	else{
		print res1 $_."\n";
	}
}
close(file);
while(my ($k,$v)=each %newhash)
{
	print res $v."\n";
}
close(res);




