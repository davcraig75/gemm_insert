#!/usr/bin/perl 
use MongoDB;
$|=1;
$VERSION="0.0200A";
$timeout=50000000;
        open(CONF,"./conf.txt") or die "Can't fiend ./conf\n";
        my %conf=();
        while (<CONF>) {
          if (/^#/) {next}
          chomp;
          if (/(.*)=(.*)/) { $conf{$1}=$2 }
          print "Conf: $1 $conf{$1} $2\n";
        }
        close (CONF);

$CONN = MongoDB::MongoClient->new(host => "mongodb://$conf{'mongodb_host'}",query_timeout =>$timeout,wtimeout =>$timeout);
$DBGENOMES = $CONN->get_database('variants');
collection_count('germline');
#collection_count('somatic');

sub collection_count {
  $q=$_[0];
  $SNVS = $DBGENOMES->get_collection($q);
  $cursor=$DBGENOMES->get_collection($q)->find();
  $c=0;

  while (my $doc=$cursor->next) {
    $chr=$doc->{'chr'};
    $hg19pos=$doc->{'hg19pos'};
    $ref=$doc->{'ref'};
    $alt=$doc->{'alt'};
    $id=$doc->{'_id'};
    if (run_cur($chr,$hg19pos,$ref,$alt,$id,$DBGENOMES,$SNVS)) {
      if ($c %100 ==0 ) { print "c:$c\n";}
      ++$c;
    }
  }
}
sub run_cur {
my ($chr,$hg19pos,$ref,$alt,$id,$DBGENOMES,$SNVS)=@_;
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

   my $cmd= Tie::IxHash->new("mapreduce"=>$q,"map" => $map, "reduce" => $reduce, 'query'=> {'chr'=>$chr,'hg19pos'=> $hg19pos,'alt'=>$alt,'ref'=>$ref},'out'=>"foo");
    my $info = $DBGENOMES->run_command($cmd);
    my $tcur=$DBGENOMES->get_collection('foo')->find_one();
    my $count=$tcur->{'value'};
    $SNVS->update({"_id"=>$id},{'$set' => {'total_count'=>int($count)}});
}
