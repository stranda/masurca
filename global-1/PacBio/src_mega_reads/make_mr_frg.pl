#!/usr/bin/env perl
#######################################
#Copyright University of Maryland 2015#
#######################################
#!/usr/bin/env perl

my $libId=$ARGV[0];
my $min_len_output=64;
$min_len_output=$ARGV[1] if($#ARGV>0);
$max_len_output=65535;

print STDOUT "{VER\n";
print STDOUT "ver:2\n";
print STDOUT "}\n";
print STDOUT "{LIB\n";
print STDOUT "act:A\n";
print STDOUT "acc:$libId\n";
print STDOUT "ori:I\n";
print STDOUT "mea:3000\n";
print STDOUT "std:300\n";
print STDOUT "src:\n";
print STDOUT ".\n";
print STDOUT "nft:3\n";
print STDOUT "fea:\n";
print STDOUT "doTrim_initialNone=1\n";
print STDOUT "doRemoveChimericReads=1\n";
print STDOUT "doRemoveSpurReads=1\n";
print STDOUT ".\n";
print STDOUT "}\n";

while($line1=<STDIN>)
{
  chomp($line1);
  if($line1 =~ /^>/)
  {
    $header=substr($line1,1);
    @f=split(/\s+/,$header);
    $readname1=substr($f[0],0,100);
    $line1=<STDIN>;
    chomp($line1);
    $len=length($line1);
    my $offset=0;
    while($offset<$len-$min_len_output){
      my $outlen = $len - $offset;
      $outlen=$max_len_output if($outlen>$max_len_output);
      $sequence1=substr($line1,$offset,$outlen);
      $clr1=0;
      $clr2=$outlen;
      print STDOUT "{FRG\n";
      print STDOUT "act:A\n";
      print STDOUT "acc:$readname1.$offset\n";
      print STDOUT "rnd:1\n";
      print STDOUT "sta:G\n";
      print STDOUT "lib:$libId\n";
      print STDOUT "pla:0\n";
      print STDOUT "loc:0\n";
      print STDOUT "src:\n.\n";
      $sequence2=$sequence1;
      $sequence2 =~ tr/Nn/Aa/;# replace Ns with As
      print STDOUT "seq:\n$sequence2\n.\n";
      $sequence1 =~ tr/ACGTNacgtn/XXXXDLLLLD/;# create fake quality scores
      print STDOUT "qlt:\n$sequence1\n.\n";
      print STDOUT "hps:\n.\n";
      print STDOUT "clv:$clr1,$clr2\n";
      print STDOUT "clr:$clr1,$clr2\n";
      print STDOUT "}\n";
      $offset+=($max_len_output-10000);
    }
  }
}
