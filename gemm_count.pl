#!/usr/bin/perl 
use MongoDB;
$|=1;
$VERSION="0.0200A";
$sample="NONE";
for ($i=0;$i<=$#ARGV;++$i) {
  if ($ARGV[$i] eq "-s") {
    ++$i;
    $sample=$ARGV[$i];
  }
}
$timeout=50000000;
        open(CONF,"/ngd-data/GEM/gemm.conf.txt") or die "Can't fiend ./conf\n";

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
collection_count('somatic');

sub collection_count {
  $q=$_[0];
  $SNVS = $DBGENOMES->get_collection($q);
  if ($sample eq "NONE") {
    $cursor=$DBGENOMES->get_collection($q)->find();
    $countsDB=$DBGENOMES->get_collection($q)->count();
  } else {
    $cursor=$DBGENOMES->get_collection($q)->find({'sample'=>$sample});
    $countsDB=$DBGENOMES->get_collection($q)->count({'sample'=>$sample});
  }
  print "sample: $sample counts: $countsDB\n";
  $c=0;
  while (my $doc=$cursor->next) {
    $chr=$doc->{'chr'};
    $hg19pos=$doc->{'hg19pos'};
    $ref=$doc->{'ref'};
    $alt=$doc->{'alt'};
    $id=$doc->{'_id'};
    $count=$DBGENOMES->get_collection($q)->count({'chr'=>$chr,'hg19pos'=> $hg19pos,'alt'=>$alt,'ref'=>$ref});
    $SNVS->update({"_id"=>$id},{'$set' => {'total_count'=>int($count)}});
    #if (run_cur($chr,$hg19pos,$ref,$alt,$id,$DBGENOMES,$SNVS)) { if ($c %1000 ==0 ) { print "c:$q $c\n";} ++$c; }
  }
}
