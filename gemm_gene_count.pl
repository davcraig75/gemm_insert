#!/usr/bin/perl 
use MongoDB;
$|=1;
$VERSION="0.0200A";
        open(CONF,"/ngd-data/GEM/gemm.conf.txt") or die "Can't fiend ./conf\n";
        my %conf=();
        while (<CONF>) {
          if (/^#/) {next}
          chomp;
          if (/(.*)=(.*)/) { $conf{$1}=$2 }
          print "Conf: $1 $conf{$1} $2\n";
        }
        close (CONF);

$CONN = MongoDB::MongoClient->new(host => "mongodb://$conf{'mongodb_host'}");
$sampleA="NONE";
for ($i=0;$i<=$#ARGV;++$i) {
  if ($ARGV[$i] eq "-s") {
    ++$i;
    $sampleA=$ARGV[$i];
    chomp($sampleA);
    print "sample: $sampleA\n";
  }
}


$DBGENOMES = $CONN->get_database('variants');
collection_count('germline');
phase('germline');
sub phase {
  $q=$_[0];
  $SNVS = $DBGENOMES->get_collection($q);
  if ($sampleA eq "NONE") {
    $cursor=$DBGENOMES->get_collection($q)->find({'project'=>'crdc',{'total_count'=>{$lt=>3},'gene_count'=>{$lt=>4}}});
    $countDB=$DBGENOMES->get_collection($q)->count({'project'=>'crdc',{'total_count'=>{$lt=>3},'gene_count'=>{$lt=>4}}});
  } else {
    $cursor=$DBGENOMES->get_collection($q)->find({'sample'=>$sampleA,{'total_count'=>{$lt=>3},'gene_count'=>{$lt=>4}}});
    $countDB=$DBGENOMES->get_collection($q)->count({'sample'=>$sampleA,{'total_count'=>{$lt=>3},'gene_count'=>{$lt=>4}}});
  }
  print "starting sample: $sampleA $countDB\n";
  $c=0;

  LOOP: while (my $doc=$cursor->next) {
  ++$c; if ($c %100 ==0 ) { print "p:$c\n";} ++$c;
    $gene_id=$doc->{'_id'};
    $SNVS->update({"_id"=>$gene_id},{'$set' => {'inheritance'=>'NA'}});
    $gene_sample=$doc->{'sample'};
    #Check eligible
    if (!($doc->{'gt1'} eq "0/0") && exists($doc->{'gt3'}) && !(exists($doc->{'gt4'}))) {
           $gene=$doc->{'gene'};
           $gene_sample=$doc->{'sample'};
           $gene_id=$doc->{'_id'};
           #Check de novo
           if (($doc->{'gt1'} eq "0/1" || $doc->{'gt1'} eq "1/1") && $doc->{'gt2'} eq "0/0" && $doc->{'gt3'} eq "0/0" && $doc->{'filter'} eq 'PASS' && $doc->{'gene_count'}<=5) {
               $SNVS->update({"_id"=>$gene_id},{'$set' => {'inheritance'=>'de_novo'}});
               next LOOP;
           }
           if (($doc->{'gt1'} eq "0/1" || $doc->{'gt1'} eq "1/1") && $doc->{'gt2'} eq "0/0" && $doc->{'gt3'} eq "0/0" && $doc->{'gene_count'}<=5) {
                $gene=$doc->{'gene'};
                $gene_sample=$doc->{'sample'};
                $gene_id=$doc->{'_id'};
               $SNVS->update({"_id"=>$gene_id},{'$set' => {'inheritance'=>'de_novo_lowQC'}});
               next LOOP;
           }
      #Check Phase
     if ($doc->{'gt1'} eq "1/1" && $doc->{'gt2'} eq "0/1" && $doc->{'gt3'} eq "0/1" && $doc->{'filter'} eq 'PASS' && $doc->{'gene_count'}<=2) {
               $gene_id=$doc->{'_id'};
               $SNVS->update({"_id"=>$gene_id},{'$set' => {'inheritance'=>'Autozygous'}});
               next LOOP;
     }
     if ($doc->{'gt1'} eq "1/1" && $doc->{'gt2'} eq "0/1" && $doc->{'gt3'} eq "0/1" && $doc->{'filter'} && $doc->{'gene_count'}<=2) {
               $gene_id=$doc->{'_id'};
               $SNVS->update({"_id"=>$gene_id},{'$set' => {'inheritance'=>'Autozygous_LowQC'}});
               next LOOP;
     }

     if ($doc->{'gt1'} eq "1/1" && ($doc->{'gt2'} eq "0/1" || $doc->{'gt3'} eq "0/1") && $doc->{'chr'} eq "chrX" && $doc->{'filter'} eq 'PASS' && $doc->{'gene_count'}<=2) {
               $gene_id=$doc->{'_id'};
               $SNVS->update({"_id"=>$gene_id},{'$set' => {'inheritance'=>'X-Linked Male'}});
               next LOOP;
     }
     if ($doc->{'gt1'} eq "1/1" && ($doc->{'gt2'} eq "0/1" || $doc->{'gt3'} eq "0/1") && $doc->{'chr'} eq "chrX" && $doc->{'filter'} eq 'PASS' && $doc->{'gene_count'}<=2) {
               $gene_id=$doc->{'_id'};
               $SNVS->update({"_id"=>$gene_id},{'$set' => {'inheritance'=>'X-Linked Male LowQC'}});
               next LOOP;
     }
      if ($doc->{'gene_count'}>=2  && $doc->{'gene_count'}<=5 ) {
        if ($doc->{'gt2'} eq "0/0"  || $doc->{'gt3'} eq "0/0") {
           $pt1=$DBGENOMES->get_collection($q)->count({'sample'=>$gene_sample,'gene'=>$gene,'gt1'=>"0/1",'gt2'=>"0/1",'gt3'=>"0/0",'filter'=>"PASS",{'total_count'=>{$lt=>3}}});
           $pt2=$DBGENOMES->get_collection($q)->count({'sample'=>$gene_sample,'gene'=>$gene,'gt1'=>"0/1",'gt2'=>"0/0",'gt3'=>"0/1",'filter'=>"PASS",{'total_count'=>{$lt=>3}}});
           if ($pt1>0 && $pt2>0) { $SNVS->update({"_id"=>$gene_id},{'$set' => {'inheritance'=>'phased_compound_het'}}); next LOOP;} 
           $pt1=$DBGENOMES->get_collection($q)->count({'sample'=>$gene_sample,'gene'=>$gene,'gt1'=>"0/1",'gt2'=>"0/1",'gt3'=>"0/0",{'total_count'=>{$lt=>3}}});
           $pt2=$DBGENOMES->get_collection($q)->count({'sample'=>$gene_sample,'gene'=>$gene,'gt1'=>"0/1",'gt2'=>"0/0",'gt3'=>"0/1",{'total_count'=>{$lt=>3}}});
           if ($pt1>0 && $pt2>0) { $SNVS->update({"_id"=>$gene_id},{'$set' => {'inheritance'=>'phased_compound_het_lowQC'}}); }
         }
       }
     }
   }
}
sub collection_count {
  $q=$_[0];
  $SNVS = $DBGENOMES->get_collection($q);
  if ($sampleA eq "NONE") {
    $cursor=$DBGENOMES->get_collection($q)->find({'project'=>'crdc',{'total_count'=>{$lt=>3}}});
    $countDB=$DBGENOMES->get_collection($q)->count({'project'=>'crdc',{'total_count'=>{$lt=>3}}});
  } else {
    $cursor=$DBGENOMES->get_collection($q)->find({'sample'=>$sampleA,{'total_count'=>{$lt=>3}}});
    $countDB=$DBGENOMES->get_collection($q)->count({'sample'=>$sampleA,{'total_count'=>{$lt=>3}}});
  }
  print "starting sample: $sampleA $countDB\n";
  $c=0;

  while (my $doc=$cursor->next) {
    $gene=$doc->{'gene'};
    $sample=$doc->{'sample'};
    $id=$doc->{'_id'};
    $count=$DBGENOMES->get_collection($q)->count({'gene'=>$gene,'sample'=> $sample, {'total_count'=>{$lt=>3}}});
    $SNVS->update({"_id"=>$id},{'$set' => {'gene_count'=>int($count)}});
    if ($c %100 ==0 ) { print "c:$c\n";} ++$c;

#    if (run_cur($gene,$sample,$id,$DBGENOMES,$SNVS)) { if ($c %1000 ==0 ) { print "c:$c\n";} ++$c; }
  }
}
sub run_cur {
my ($gene,$sample,$id,$DBGENOMES,$SNVS)=@_;
my $map = <<MAP;
  function() {
       emit(this.id,1);
  }
MAP

my $reduce = <<REDUCE;
    function(prev,current) {
        count = 0;
        current.forEach(function(item) {
            count = count+1;
        });
        return count;
    }
REDUCE

    my $cmd= Tie::IxHash->new("mapreduce"=>$q,"map" => $map, "reduce" => $reduce, 'query'=> {'gene'=>$gene,'sample'=> $sample, {'total_count'=>{$lt=>3}}},'out'=>"foo");
    my $info = $DBGENOMES->run_command($cmd);
    my $tcur=$DBGENOMES->get_collection('foo')->find_one();
    my $count=$tcur->{'value'};
    $SNVS->update({"_id"=>$id},{'$set' => {'gene_count'=>int($count)}});
}
