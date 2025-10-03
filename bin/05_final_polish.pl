#!/usr/bin/env perl
$INPUT=$ARGV[0];    # final_sr2lr_sv_vcf
$INPUT2=$ARGV[1];   # shared.polished.bedpe
$OUTPUT=$ARGV[2];   # final.sr2lr.sv.bedpe 
$OUTPUT2=$ARGV[3];  # final.sr2lr.polished.vcf

#$shared=$INPUT;
#@filelist=split(/\//,$INPUT);
#$filelist[-1]="shared.polished.bedpe";
#$newpath=join("/",@filelist);
open(file,"$INPUT2");

while(<file>)
{
	chomp;
	@token=split(/\t/,$_);
	$lr2sr{$token[10]}=0;
	# print $token[10]."\n";
}
close(file);

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
		$oldid=$token[2];
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
		if($token[9] ne "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN")
		{
			if(exists($lr2sr{$oldid}))
			{
				$seqflag="LR-SR";
			}else{
				$seqflag="LR";	
			}
			$tags=$token[9];
		}
		if($token[10] ne "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN")
		{
			$seqflag="SR";
			$tags=$token[10];
		}

		if(grep(/\>/,$token[4]))
		{
			$ts=$token[4];
			$ts=~s/\>//g;
			$ts=~s/\<//g;
			$type=$ts;
		}
		$loc1_1=$loc1-1;
		$loc1_2=$loc1+1;
		$loc2_1=$loc2-1;
		$loc2_2=$loc2+1;
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
				# print res1 $token[0]."\t".$token[1]."\t".$token[2]."\t".$token[3]."\t"."<".$type.">"."\t".$token[5]."\t".$token[6]."\t".$resline."\t".$token[8]."\t".$token[9]."\t".$token[10]."\n";
				print res1 $token[0]."\t".$token[1]."\t".$token[2]."\t".$token[3]."\t"."<".$type.">"."\t".$token[5]."\t".$token[6]."\t".$resline."\t".$token[8]."\t".$tags."\n";
			}else{
				@t1=split(/SVTYPE=/,$token[7]);
				@t2=split(/\;SVME/,$t1[1]);
				$resline=$t1[0]."SVTYPE=".$type.";SVME".$t2[1];
				print res1 $token[0]."\t".$token[1]."\t".$token[2]."\t".$token[3]."\t"."<".$type.">"."\t".$token[5]."\t".$token[6]."\t".$resline."\t".$token[8]."\t".$tags."\n";
			}
		}else{
			@t1=split(/SVTYPE=/,$token[7]);
			@t2=split(/\;SVME/,$t1[1]);
			$resline=$t1[0]."SVTYPE=".$type.";SVME".$t2[1];
			print res1 $token[0]."\t".$token[1]."\t".$token[2]."\t".$token[3]."\t".$token[4]."\t".$token[5]."\t".$token[6]."\t".$resline."\t".$token[8]."\t".$tags."\n";
			# print res1 $token[0]."\t".$token[1]."\t".$token[2]."\t".$token[3]."\t"."<".$type.">"."\t".$token[5]."\t".$token[6]."\t".$resline."\t".$token[8]."\t".$token[9]."\t".$token[10]."\n";
		}

		if($type eq "INV")
		{
			$len=abs($loc2-$loc1);
		}


		if(!exists($hash{$chr1."\t".$loc1}) and !exists($hash{$chr1."\t".$loc1_1}) and !exists($hash{$chr1."\t".$loc1_2}) and !exists($hash{$chr2."\t".$loc2}) and !exists($hash{$chr2."\t".$loc2_1}) and !exists($hash{$chr2."\t".$loc2_2}))
		{
			@token_1=split(/STRANDS=/,$_);
			@token_2=split(/\;/,$token_1[1]);
			@token_3=split(//,$token_2[0]);
			$strand1=$token_3[0];
			$strand2=$token_3[1];
			$hash{$chr1."\t".$loc1}=$token[2];
			$hash{$chr2."\t".$loc2}=$token[2];
			$hash{$chr1."\t".$loc1_1}=$token[2];
			$hash{$chr1."\t".$loc1_2}=$token[2];
			$hash{$chr2."\t".$loc2_1}=$token[2];
			$hash{$chr2."\t".$loc2_1}=$token[2];
			$newhash{$token[2]}=$chr1."\t".$loc1."\t".$loc1."\t".$chr2."\t".$loc2."\t".$loc2."\t".$type."\t".$id."\t".$len."\t".$seqflag."\t".$strand1."\t".$strand2;

		}else{
			if($type eq "INV")
			{
				if(exists($hash{$chr1."\t".$loc1}) and  grep(/SUPP=2/,$_))
				{
					$newk=$hash{$chr1."\t".$loc1};
					delete $newhash{$newk};
					$newhash{$token[2]}=$chr1."\t".$loc1."\t".$loc1."\t".$chr2."\t".$loc2."\t".$loc2."\t".$type."\t".$id."\t".$len."\t".$seqflag."\t".$strand1."\t".$strand2;
				}
				if(exists($hash{$chr1."\t".$loc1_1}) and  grep(/SUPP=2/,$_))
				{
					$newk=$hash{$chr1."\t".$loc1_1};
					delete $newhash{$newk};
					$newhash{$token[2]}=$chr1."\t".$loc1."\t".$loc1."\t".$chr2."\t".$loc2."\t".$loc2."\t".$type."\t".$id."\t".$len."\t".$seqflag."\t".$strand1."\t".$strand2;
				}
				if(exists($hash{$chr1."\t".$loc1_2}) and  grep(/SUPP=2/,$_))
				{
					$newk=$hash{$chr1."\t".$loc1_2};
					delete $newhash{$newk};
					$newhash{$token[2]}=$chr1."\t".$loc1."\t".$loc1."\t".$chr2."\t".$loc2."\t".$loc2."\t".$type."\t".$id."\t".$len."\t".$seqflag."\t".$strand1."\t".$strand2;
				}
				if(exists($hash{$chr2."\t".$chr2}) and  grep(/SUPP=2/,$_))
				{
					$newk=$hash{$chr2."\t".$chr2};
					delete $newhash{$newk};
					$newhash{$token[2]}=$chr1."\t".$loc1."\t".$loc1."\t".$chr2."\t".$loc2."\t".$loc2."\t".$type."\t".$id."\t".$len."\t".$seqflag."\t".$strand1."\t".$strand2;
				}
				if(exists($hash{$chr2."\t".$loc2_1}) and grep(/SUPP=2/,$_))
				{
					$newk=$hash{$chr2."\t".$loc2_1};
					delete $newhash{$newk};
					$newhash{$token[2]}=$chr1."\t".$loc1."\t".$loc1."\t".$chr2."\t".$loc2."\t".$loc2."\t".$type."\t".$id."\t".$len."\t".$seqflag."\t".$strand1."\t".$strand2;
				}
				if(exists($hash{$chr2."\t".$loc2_2}) and grep(/SUPP=2/,$_))
				{
					$newk=$hash{$chr2."\t".$loc2_2};
					delete $newhash{$newk};
					$newhash{$token[2]}=$chr1."\t".$loc1."\t".$loc1."\t".$chr2."\t".$loc2."\t".$loc2."\t".$type."\t".$id."\t".$len."\t".$seqflag."\t".$strand1."\t".$strand2;
				}
			}
		}


		# }
	}
	else{
		if(grep(/CHROM/,$_))
		{
			$out1=$_;
			$out1=~s/	NULL//g;
			print res1 $out1."\n";
		}else{
			print res1 $_."\n";
		}
		
	}
}
close(file);
while(my ($k,$v)=each %newhash)
{
	print res $v."\t".$k."\n";
}
close(res);
