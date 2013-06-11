#!/usr/bin/perl -w 
use threads;
use threads::shared;
use MongoDB;
$MAX_THREADS=1;
$|=1;
$VERSION="0.0200A";
$host="mongodb://spock.tgen.org:27017";
$CONN = MongoDB::MongoClient->new(host => $host);
$DBGENOMES = $CONN->get_database('variants');
collection_count('germline');
collection_count('somatic');

sub collection_count {
  $q=$_[0];
  $SNVS = $DBGENOMES->get_collection($q);
  $cursor=$DBGENOMES->get_collection($q);
  $c=0;

  while (my $doc=$cursor->next) {
    $gene=$doc->{'gene'};
    $sample=$doc->{'sample'};
    $id=$doc->{'_id'};
    if (run_cur($gene,$sample,$id,$DBGENOMES,$SNVS)) {
      if ($c %10000 ==0 ) { print "c:$c\n";}
      ++$c;
    }
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

   my $cmd= Tie::IxHash->new("mapreduce"=>$q,"map" => $map, "reduce" => $reduce, 'query'=> {'gene'=>$gene,'sample'=> $sample, 'count'=>int(1)},'out'=>"foo");
    my $info = $DBGENOMES->run_command($cmd);
    my $tcur=$DBGENOMES->get_collection('foo')->find_one();
    my $count=$tcur->{'value'};
    print "$gene $sample $count\n";
    $SNVS->update({"_id"=>$id},{'$set' => {'gene_count'=>int($count)}},{'upsert' => 1});
}
