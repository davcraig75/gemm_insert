#!/usr/bin/perl 
use strict;
use Time::HiRes;
use Time::HiRes qw(usleep);
use Storable qw(dclone);
use File::stat;
use Time::localtime;
use threads;
use threads::shared;
use MongoDB;
# Version info
## GLOBAL VARIABLES
	$MongoDB::BSON::use_binary = 4;
	$|=1;
        our $VERSION="0.07000";
	our $GEM_PATH="/ngd-data/Craig_lab/dcraig/gedi";
	open(CONF,"./conf.txt") or die "Can't fiend ./conf\n";
        my %conf=();
        while (<CONF>) {
          if (/^#/) {next}
          chomp;
          if (/(.*)=(.*)/) { $conf{$1}=$2 }
          print "Conf: $1 $conf{$1} $2\n";
        }
        close (CONF);
        my $ANNOTATION_SERVER=$conf{'annotation_host'};
        my $SERVER=$conf{'mongodb_host'};
	open (OUT,">$GEM_PATH/gemm_insert.$SERVER.out") or die "cant open log file";
	print OUT "------gemm_insert version: $VERSION server: $SERVER-------\n";
	our $MAX_THREADS=3;
	my $timeout=10000000; 
	our $CONN = MongoDB::Connection->new(host => "$SERVER");
	our $CONNA = MongoDB::Connection->new(host => "$ANNOTATION_SERVER");
	$CONN->query_timeout($timeout); $CONN->wtimeout($timeout);
	our $DBGENOMES = $CONN->get_database('variants');
	our $DBSAMPLES = $CONN->get_database('samples');
	our $ANNOTATE = $CONNA->get_database('annotation');
	our $MONGO_FILES = $DBGENOMES->get_collection('files');
	our $PATIENTS = $DBSAMPLES->get_collection('patients');
	our $GENES = $ANNOTATE->get_collection('genes');
	our $LEVEL2G = $ANNOTATE->get_collection('level2');
	our $DBSNP = $ANNOTATE->get_collection('dbsnp');
	our $DBNSFP = $ANNOTATE->get_collection('dbNSFP');
	our $CLINVARDB = $ANNOTATE->get_collection('clinvar');
	our $COSMICDB = $ANNOTATE->get_collection('cosmic');
	our $TCGADB = $ANNOTATE->get_collection('tcga');
	our $CHECK = $DBGENOMES->get_collection('somatic');

## MAIN  

	open (VCF_PATHS,"$GEM_PATH/$conf{'vcf_path'}") or die "Can't open $GEM_PATH/$conf{'vcf_path'}\n";
	while (<VCF_PATHS>) {
		chomp;
		my $vcf_path=$_;
		print OUT "-VCF Folder path: $vcf_path\n";
		my @list=`find $vcf_path -name "*.vcf"`;
		&clean_files(\@list,$MONGO_FILES, $DBGENOMES);
		while (my $file=shift(@list)) {
			chomp($file); 
			if(threads->new(\&process_folder,$file)==1) {print "+File processing: $file\n"; } #else {print "-File skipped: $file\n";}
#			if(process_folder($file)) {print "+File processing: $file\n"; } #else {print "-File skipped: $file\n";}
			while (threads->list(threads::running) >= $MAX_THREADS) { sleep 10; }
		}
	}
	while (threads->list(threads::running) > 0) { sleep 10; }
	foreach my $thr (threads->list(threads::joinable)) {
		$thr->join;
	}
	close (VCF_PATHS);
	print "Done.\n";
	close (OUT);
        my $SNVS = $DBGENOMES->get_collection('somatic');
        indexing($SNVS);
        my $SNVS = $DBGENOMES->get_collection('germline');
        indexing($SNVS);



## PROCESS FILE
sub process_folder {
    my $file=$_[0];
    my %file_info=();
    my @threads_ar=();
    $file_info{'filepath'}=$file;
    if ($file_info{'filepath'}=~/.*VCF\/(.*?)\/(.*?)\/(.*?)\/(.*.vcf)$/) {
      $file_info{'project'}=$1;
      $file_info{'sample'}="$1.$2";
      $file_info{'tissue'}=$3;
      $file_info{'filename'}=$4;
    } elsif ($file_info{'filepath'}=~/.*VCF\/(.*?)\/(.*?)\/(.*.vcf)$/) {
      $file_info{'project'}=$1;
      $file_info{'sample'}="$1.$2";
      $file_info{'tissue'}="germline";
      $file_info{'filename'}=$3;
    } elsif ($file_info{'filepath'}=~/.*VCF\/(.*?)\/(.*.vcf)$/) {
      $file_info{'project'}=$1;
      $file_info{'sample'}="$1.$2";
      $file_info{'tissue'}="germline";
      $file_info{'filename'}=$2;
    } elsif ($file_info{'filepath'}=~/.*VCF\/(.*.vcf)$/) {
      $file_info{'project'}="no_project";
      $file_info{'sample'}=$1;
      $file_info{'tissue'}="germline";
      $file_info{'filename'}=$1;
    }
    my $SNVS = $DBGENOMES->get_collection('somatic');
    $file_info{'database'}="somatic";
    if ($file_info{'tissue'}=~/germ/i) {
      $file_info{'tissue'}="germline";
      $file_info{'database'}="germline";
      $SNVS = $DBGENOMES->get_collection('germline');
      #duggan test
      #my $temp=`head -100 $file | grep "#CHROM"`;
      #chomp ($temp);
      #my @temp_ar=split(/\t/,$temp);
      #$file_info{'sample'}=$temp_ar[9];
    } 
    #Timestamp
    my $date = ctime(stat($file)->mtime);
    $file_info{'date'}=$date;
    my $cur_file=$MONGO_FILES->find({'filename'=>$file_info{'filename'}});
    if (my $col=$cur_file->next) { 
       if ($col->{'date'} ne $file_info{'date'}) {
         $MONGO_FILES -> remove({'filename'=>$file_info{'filename'}});
         $SNVS -> remove({'filename'=>$file_info{'filename'}});
       } else {
         return 0;
       }
    }
    $SNVS -> remove({'filename'=>$file_info{'filename'}});
    #run_snpeff($file_info{'filepath'});
    #threads->new(\&parse_vcf,\%file_info,$SNVS,$MONGO_FILES);
    print "-preparing $file_info{'project'},$file_info{'sample'},$file_info{'tissue'},$file_info{'filename'}\n";
    parse_vcf(\%file_info,$SNVS,$MONGO_FILES);
    $MONGO_FILES -> insert(\%file_info);
    #`mkdir -p $GEM_PATH/backup/VCF/$file_info{'project'}/$file_info{'sample'}/$file_info{'tissue'}`;
    #`cp -f $file $GEM_PATH/backup/VCF/$file_info{'project'}/$file_info{'sample'}/$file_info{'tissue'}/$file_info{'filename'}&`;
	return 1;
}


##################################
# process files
##################################
sub parse_vcf {
	my $file_info=$_[0]; 
	my $SNVS=$_[1];
	my $MONGO_FILES=$_[2];
        #my $idx = Tie::IxHash->new('filename');
        #$SNVS->ensure_index($idx);
	open (FILE, "$file_info->{'filepath'}") or die "Can't open file: $file_info->{'filepath'}!\n";
	my ($line,$i,$j,$set3);
	my %repeat=();
        my $counter=0;
	HEADER: while ($line=<FILE>) {
		chomp($line);
		if (!($line=~/^\#/)) {last HEADER}
	}
	MAIN_LOOP: while ($line) {
		chomp($line);
                # Added due to fusions
                $line=~s/FUSED_GENE=/GENE=/g;
                $line=~s/GENE1=/GENE=/g;
                $line=~s/GENE2=/GENE=/g;
                ++$counter;
                #if ($counter%100000 == 0) {print "-Counter $file_info->{'filepath'} $counter\n";}
		## Parsing Mandatory part, and genotypes
		$line=~s/\"//g;
		my @temp=split(/\t+/,$line);
		$line=(<FILE>); 
                #if ($line=~/EFF/) {}else{next MAIN_LOOP}
		my %var_info= (
		'chr'=>$temp[0],
		'hg19pos'=>int($temp[1]),
		'aberration_name'=>$temp[2],
		'ref'=>$temp[3],
		'alt'=>$temp[4],
		'qual'=>$temp[5],
		'aberration_filter'=>$temp[6],
		'filter'=>$temp[6],
		);
    my $var_info_ref=\%var_info;
    my $info=$temp[7];
    if ($var_info{'chr'}!~/chr/) {$var_info{'chr'}=join("","chr",$var_info{'chr'})}
    if ($var_info{'qual'} eq ".") {$var_info{"qual"} =255;$var_info{'aberration_filter'}='PASS';}
    #if ($var_info{'aberration_filter'} eq '.') {$var_info{'aberration_filter'}='PASS';}
    if ($var_info{'aberration_filter'} ne 'PASS') {next MAIN_LOOP}
    if($var_info{'qual'}>3) {$var_info{'qual'}=int($var_info{'qual'})};
	$var_info{'coordinate'}=join(":",$var_info{'chr'},$var_info{'hg19pos'},$var_info{'ref'},$var_info{'alt'});
    $var_info{'deg'}=0;
    chomp($info);
    # Determine and parse multiple genotypes for germline
    if (exists($temp[9])) { 
      my @gt=@temp[9..$#temp]; 
      my @fields=split(/:/,lc($temp[8]));
      for ($i=0;$i<=$#gt;++$i) {
        my @temp4=split(/:/,$gt[$i]);
        for ($j=0;$j<=$#temp4;++$j) {
           $var_info{join("",$fields[$j],$i+1)}=$temp4[$j];
        }
      }
    }
    $var_info_ref=add_in_hash($var_info_ref,$file_info);
    if (parse_info($var_info_ref,$info) == -1 ) {next MAIN_LOOP}
    if ($var_info_ref->{'alt'} eq "<SV>") { if ($var_info_ref->{'sv_gap'}<2200) { next MAIN_LOOP} }
    if (assign_aberration($var_info_ref) == -1) {next MAIN_LOOP}
    foreach my $t (@{$var_info_ref->{'genes'}}) {
      my $anom_var_info_ref=dclone $var_info_ref;
      $anom_var_info_ref->{'gene'}=$t;
      if (annotate_line($anom_var_info_ref) == -1 ) {next MAIN_LOOP}
      if ($anom_var_info_ref->{'aberration_filter'} eq "PASS") {
        if (insert_document($anom_var_info_ref,\%repeat,$SNVS)== -1) {next MAIN_LOOP}
      }
    }
  }
  close (FILE);
  print OUT "Added project: $file_info->{'project'}  sample: $file_info->{'sample'}  tissue: $file_info->{'tissue'}  filename:$file_info->{'filename'}\n";
}
##################
###   assign_aberration
######################

sub assign_aberration {
  my $var_info_ref=$_[0];
  if (ref($var_info_ref) eq "HASH") {
    if ($var_info_ref->{'alt'} eq "<SV>") {
      $var_info_ref->{'aberration_type'}="STRUCTURAL_VARIANT";
      if ($var_info_ref->{'qual'}<6) {
        $var_info_ref->{'aberration_filter'}="LOW";
      } else {
        $var_info_ref->{'aberration_filter'}="PASS";
      }
      $var_info_ref->{'aberration_value'}=$var_info_ref->{'qual'};
      $var_info_ref->{'aberration_effect'}="DELETERIOUS";
      $var_info_ref->{'biomarker_type'}="gene";
    } elsif ($var_info_ref->{'alt'} eq "<DEL>" || $var_info_ref->{'alt'} eq "<LOSS>") {
      $var_info_ref->{'aberration_type'}="CNV_LOSS";
      $var_info_ref->{'aberration_value'}=$var_info_ref->{'log2fc'};
      $var_info_ref->{'aberration_effect'}="DELETION";
      $var_info_ref->{'aberration_name'}=$var_info_ref->{'chr'} . ":" . $var_info_ref->{'hg19pos'} . "-" . $var_info_ref->{'end'} ;
      $var_info_ref->{'biomarker_type'}="gene";
      if ($var_info_ref->{'qual'}<1) {$var_info_ref->{'aberration_filter'}='FAIL'}
    } elsif ($var_info_ref->{'alt'} eq "<DUP>" || $var_info_ref->{'alt'} eq "<GAIN>") {
      $var_info_ref->{'aberration_type'}="CNV_GAIN";
      $var_info_ref->{'aberration_value'}=$var_info_ref->{'log2fc'};
      $var_info_ref->{'aberration_effect'}="AMPLIFICATION";
      $var_info_ref->{'aberration_name'}=$var_info_ref->{'chr'} . ":" . $var_info_ref->{'hg19pos'} . "-" . $var_info_ref->{'end'} ;
      $var_info_ref->{'biomarker_type'}="gene";
      if ($var_info_ref->{'qual'}<1) {$var_info_ref->{'aberration_filter'}='FAIL'}
    } elsif ($var_info_ref->{'alt'} eq "<FUSION>") {
      $var_info_ref->{'aberration_type'}="FUSED_GENE";
      $var_info_ref->{'aberration_name'}=$var_info_ref->{'fused_gene'};
      $var_info_ref->{'gene'}=$var_info_ref->{'fused_gene'};
      $var_info_ref->{'aberration_effect'}="FUSION";
      $var_info_ref->{'aberration_value'}=$var_info_ref->{'qual'};
      $var_info_ref->{'aberration_filter'}="PASS";
    } elsif ($var_info_ref->{'alt'} eq "<EXP>") {
      $var_info_ref->{'aberration_value'}=sprintf("%1.2f",$var_info_ref->{'log2fc'});
      if ($var_info_ref->{'aberration_value'}>0) {
         $var_info_ref->{'aberration_type'}="EXP_OVER";
         $var_info_ref->{'aberration_effect'}="OVER_EXPRESSED";
         $var_info_ref->{'biomarker_type'}="gene";
      }elsif ($var_info_ref->{'aberration_value'}<0) {
         $var_info_ref->{'aberration_effect'}="UNDER_EXPRESSED";
         $var_info_ref->{'aberration_type'}="EXP_UNDER";
         $var_info_ref->{'biomarker_type'}="gene";
      }
    } elsif (length($var_info_ref->{'ref'}) >1) {
      $var_info_ref->{'aberration_type'}="SMALL_DELETION";
      $var_info_ref->{'aberration_effect'}=$var_info_ref->{'functional_change'};
      $var_info_ref->{'biomarker_type'}="gene";
      $var_info_ref->{'aberration_name'}=$var_info_ref->{'functional_change'};
      if ($var_info_ref->{'aberration_effect'}=~/(CODON)|(FRAME)|(SPLICE)|(START_)|(STOP_)|(NON_SYNONYMOUS)/) { $var_info_ref->{'aberration_filter'}='PASS'; } else { $var_info_ref->{'aberration_filter'}='NON_GENIC'; }
    } elsif (length($var_info_ref->{'alt'})>1 ) {
      $var_info_ref->{'aberration_type'}="SMALL_INSERTION";
      $var_info_ref->{'biomarker_type'}="gene";
      $var_info_ref->{'aberration_effect'}=$var_info_ref->{'functional_change'};
      $var_info_ref->{'aberration_name'}=$var_info_ref->{'functional_change'};
      if ($var_info_ref->{'aberration_effect'}=~/(CODON)|(FRAME)|(SPLICE)|(START_)|(STOP_)|(NON_SYNONYMOUS)/) { $var_info_ref->{'aberration_filter'}='PASS'; } else { $var_info_ref->{'aberration_filter'}='NON_GENIC'; }
    } elsif ($var_info_ref->{'type'} eq "SNV") {
      $var_info_ref->{'aberration_type'}="SNV";
      $var_info_ref->{'biomarker_type'}="gene";
      $var_info_ref->{'aberration_effect'}=$var_info_ref->{'functional_change'};
      if ($var_info_ref->{'aberration_effect'}=~/(CODON)|(FRAME)|(SPLICE)|(START_)|(STOP_)|(NON_SYNONYMOUS)/) { $var_info_ref->{'aberration_filter'}='PASS'; } else { $var_info_ref->{'aberration_filter'}='NON_GENIC'; }
    }
    if ($var_info_ref->{'aberration_type'} eq 'SNV' || $var_info_ref->{'aberration_type'} eq 'SMALL_INSERTION' || $var_info_ref->{'aberration_type'} eq 'SMALL_DELETION') {
      if (exists($var_info_ref->{'amino_acid_change'})) {
        $var_info_ref->{'aberration_value'}=$var_info_ref->{'amino_acid_change'};
      } else {
        $var_info_ref->{'aberration_value'}=$var_info_ref->{'functional_change'};
      }
    }
    if ($var_info_ref->{'aberration_name'} eq ".") { $var_info_ref->{'aberration_name'}=$var_info_ref->{'chr'} . ":" . $var_info_ref->{'hg19pos'}; }
  }
  return 1;
}

##################
##### check_snpeff
##################
sub run_snpeff{
  my $file=$_[0];
  my $temp=`head -500 $file > ./temp`;
  my @grep=`grep "SnpEff 3" temp`;
  return 1; 
  if (scalar(@grep)>0) {return 1}
  @grep=`grep UnifiedGenotyper temp`;
  if (scalar(@grep)>0) {return 1}
  @grep=`grep omatic temp`;
  if (scalar(@grep)>0) {
    print "-snpeff converting $file\n";
    `rm -f temp`;
    `cp $file temp`; 
    `java -jar ./misc/snpEff.jar eff -no-downstream  -no-intergenic -no-intron -no-upstream -minQ 20 -canon GRCh37.64 -c ./misc/snpEff.config -o vcf ./temp >  $file`;
  }
  return 1; 
}

##################
####   PARSE INFO
#######################
sub parse_info {
  my $var_info_ref=$_[0];
  my $info=$_[1];
  $var_info_ref->{"genes"}=[];
  #bless($var_info_ref->{"genes"},"MongoDB::BSON::Array");
  my @words=split(/\;+/,$info);
  LOOP2: foreach my $word (@words) {
      my ($key,$val)="";
      if ($key eq "FUSED_GENE") {$key="gene"}
      if ($word=~/=/) { ($key,$val)=split(/=/,$word) }else { next LOOP2 }
      if ($key eq "PILEUP") { next LOOP2} 
      if ($key eq "TYPE") { $val=~s/somatic_//g;} 
      if ($key eq "TYPE" && $val=~/LOH/) { return -1}
      if ($key eq "EFF") {
        if ($var_info_ref->{'qual'} <20) {return -1}
        if (check_snpeff($var_info_ref,$key,$val)==-1) {return -1}
      } elsif ($key eq "GENE" || $key eq "gene") {
        $key=lc($key);
        $var_info_ref->{$key}=($val);
        push(@{$var_info_ref->{'genes'}},$val);
      } else {
        $key=lc($key);
        $var_info_ref->{$key}=($val);
      }
   }
   if (($var_info_ref->{'alt'} eq "<SV>" || $var_info_ref->{'alt'} eq "<FUSION>") && scalar(@{$var_info_ref->{'genes'}})<2) {push(@{$var_info_ref->{'genes'}},"NONE")}
   return 1;
}
##################
##   insert_document
#####################
sub insert_document {
  my $var_info_ref=$_[0];
  my $repeat_ref=$_[1];
  my $SNVS=$_[2];
  my $id="";
  my $cur="";
  if (ref($var_info_ref) eq "HASH") {
    if (exists($repeat_ref->{"$var_info_ref->{'aberration_name'}-$var_info_ref->{'gene'}" })) {
          return -1;
    } else {
        $repeat_ref->{"$var_info_ref->{'aberration_name'}-$var_info_ref->{'gene'}" }=1;
    }
    $var_info_ref->{'snpeff_insert'}=$VERSION;
    if (insert_check($var_info_ref)==1) {
      $var_info_ref->{'rare_variant'}=int(1);
      $var_info_ref->{'chr'}="$var_info_ref->{'chr'}";
      $var_info_ref->{'hg19pos'}="$var_info_ref->{'hg19pos'}";
      $var_info_ref->{'count'}=int(1);
      $var_info_ref->{'total_count'}=int(1);
      if (exists($var_info_ref->{'dp1'})){ $var_info_ref->{'dp1'}=int($var_info_ref->{'dp1'}) }
      if (exists($var_info_ref->{'dp2'})){$var_info_ref->{'dp2'}=int($var_info_ref->{'dp2'})}
      if (exists($var_info_ref->{'dp3'})){$var_info_ref->{'dp3'}=int($var_info_ref->{'dp3'})}
      if (exists($var_info_ref->{'dp4'})){$var_info_ref->{'dp4'}=int($var_info_ref->{'dp4'})}
      if (exists($var_info_ref->{'dp5'})){$var_info_ref->{'dp5'}=int($var_info_ref->{'dp5'})}
      if (exists($var_info_ref->{'qual'})){$var_info_ref->{'qual'}=int($var_info_ref->{'qual'})}
      $var_info_ref->{'hg19pos'}="$var_info_ref->{'hg19pos'}";
      $id=$SNVS->insert($var_info_ref, {safe => 1});
    }
  }
   return 1;
}

##################
#   annotation via database 
####################
sub annotate_line {
  my $var_info_ref=$_[0];
  if ($var_info_ref->{'functional_change'}) {
    $var_info_ref=add_in_hash($var_info_ref,query_DBSNP($var_info_ref->{'chr'},$var_info_ref->{'hg19pos'},$var_info_ref->{'alt'},$var_info_ref->{'ref'}));
    $var_info_ref=add_in_hash($var_info_ref,CLINVAR($var_info_ref->{'chr'},$var_info_ref->{'hg19pos'},$var_info_ref->{'alt'},$var_info_ref->{'ref'}));
    $var_info_ref=add_in_hash($var_info_ref,COSMIC($var_info_ref->{'chr'},$var_info_ref->{'hg19pos'},$var_info_ref->{'alt'},$var_info_ref->{'ref'}));
    if ($var_info_ref->{'deg'}>2) { 
        $var_info_ref=add_in_hash($var_info_ref,NSFP($var_info_ref->{'chr'},$var_info_ref->{'hg19pos'},$var_info_ref->{'alt'},$var_info_ref->{'ref'}));
     }
    $var_info_ref->{'count'}=int(1);
   }
   if (exists($var_info_ref->{'gene'} )) { 
      $var_info_ref=add_in_hash($var_info_ref,IGENE($var_info_ref->{'gene'}));
      $var_info_ref=add_in_hash($var_info_ref,LEVEL2($var_info_ref->{'gene'}));
   }
   # su2c specific
   if ($var_info_ref->{'sample'}=~/su2c/) {$var_info_ref=add_in_hash($var_info_ref,SAMPLEIN($var_info_ref->{'sample'}))} # su2c specifici
   if ($var_info_ref->{'database'} eq "somatic") { $var_info_ref=add_in_hash($var_info_ref,TCGA($var_info_ref->{'chr'},$var_info_ref->{'hg19pos'},$var_info_ref->{'alt'},$var_info_ref->{'ref'})); }
  return 1;
}
##################
##   add in hash
#####################
sub add_in_hash {
  my $main_hash_ref=$_[0];
  my $sub_hash_ref=$_[1];
  foreach my $key (keys %$sub_hash_ref) {
    unless (exists ($main_hash_ref->{$key})) {
      $main_hash_ref->{$key}=$sub_hash_ref->{$key};
    }
  }
  return $main_hash_ref;
}
sub replace_in_hash {
  my $main_hash_ref=$_[0];
  my $sub_hash_ref=$_[1];
  foreach my $key (keys %$sub_hash_ref) {
      $main_hash_ref->{$key}=$sub_hash_ref->{$key};
  }
  return $main_hash_ref;
}

##################
# Decide whether to insert
# ##################
#
sub insert_check {
   my $hash_ref=$_[0];
   my $toinsert=0;
   if (exists($hash_ref->{'dp1'})) {
     if ($hash_ref->{'dp1'}<10) {return 0}
   }
   foreach my $key (keys %{$hash_ref}) { 
     my $v=$hash_ref->{$key};
     $hash_ref->{$key}="$v"; 
    }
   if (exists($hash_ref->{'dp1'})) { if ($hash_ref->{'dp1'}<10) { return 0} }
#   if (!(exists($hash_ref->{'1000G_ASN_AF'}))) { $hash_ref->{'1000G_ASN_AF'}=-1.01}
#   if (!(exists($hash_ref->{'ESP6500_AA_AF'}))) { $hash_ref->{'ESP6500_AA_AF'}=-1.01} 
#   if (!(exists($hash_ref->{'ESP6500_EA_AF'}))) { $hash_ref->{'ESP6500_EA_AF'}=-1.01} 
   if ($hash_ref->{'esp6500_ea_af'}>0.03 || $hash_ref->{'esp6500_aa_af'}>0.03 ||  $hash_ref->{'1000gp1_af'}>0.03 ) { $toinsert=0; return $toinsert }

   MongoDB::force_double($hash_ref->{'1000gp1_af'});
   MongoDB::force_double($hash_ref->{'esp6500_ea_af'});
   MongoDB::force_double($hash_ref->{'esp6500_aa_af'});
   if ($hash_ref->{'aberration_value'} eq "FAIL"){ $toinsert=0;}
   if (($hash_ref->{'aberration_type'} eq "CNV_GAIN" ||$hash_ref->{'aberration_type'} eq "CNV_LOSS") && ($hash_ref->{'aberration_value'} <1 && $hash_ref->{'aberration_value'} >-1)){$toinsert=0;}
   if (exists ($hash_ref->{'gene'})) { $toinsert=1; }
   if ( exists ($hash_ref->{'functional_change'})) { if($hash_ref->{'functional_change'}=~/(CODON)|(FRAME)|(SPLICE)|(NON_SYNONYMOUS)|(START_)|(STOP_)|(SYNONYMOUS)|(UTR)|(EXON)/) { $toinsert=1; }else {$toinsert=0} }
   return $toinsert;
}

sub indexing {
  my $SNVS = $_[0];
	my ($idx,$val);
	my @vals=('coordinate','sample','chr','hg19pos','ESP6500_AA_AF','ESP6500_EA_AF','ref','alt','dbsnp',
	'aberration_type','aberration_name','aberration_value','aberration_effect','functional_change','CLNSRC','CLNDBN','CLNDSDB','file',
	'gene','count','total_count','GENEINFO','KGPhase1','OTHERKG','ac','af','baseranksum','clinvar','codon_change','coordinate','count','total_count',
        'project','panels');
	foreach $val (@vals) {
		$idx = Tie::IxHash->new($val => 1);
		$SNVS->ensure_index($idx);
	}
	$idx = Tie::IxHash->new("chr" => 1,"hg19pos" => 1, "alt" => 1, "ref" => -1); 
	$SNVS->ensure_index($idx);
}

sub clean_files {
 my @list=@{$_[0]};
 my $MONGO_FILES=$_[1];
 my $DBGENOMES=$_[2];
 my %filename_hash=();
 while (my $filename=shift(@list)) {
   chomp($filename); 
   $filename=~s/.*\///;
   $filename_hash{"$filename"}=1;
 }
 my $cur=$MONGO_FILES->find();
 while (my $doc=$cur->next) {
     my $mongo_filename=$doc->{'filename'};
     my $tissue=$doc->{'tissue'};
     if (!(exists($filename_hash{"$mongo_filename"}))) {
         print "Scrubbing: $mongo_filename\n";
         $MONGO_FILES -> remove({'filename'=>"$mongo_filename"});
         my $lSNVS = $DBGENOMES->get_collection('somatic');
         if ($tissue eq "germline") { my $lSNVS = $DBGENOMES->get_collection('germline'); } 
         $lSNVS -> remove({'filename'=>"$mongo_filename"});
     }
 }
}
sub NSFP {
   my $chr=$_[0];my $hg19pos=int($_[1]);my $alt=$_[2];my $ref=$_[3];
   my %temp_h=();
   $chr=~s/chr//g;
   my $cur=$DBNSFP->find({'chr'=>"$chr",'hg19pos'=>int($hg19pos),'alt'=>$alt,'ref'=>$ref});
   if (my $col=$cur->next) {%temp_h=%$col; delete $temp_h{"chr"};delete $temp_h{"_id"};delete $temp_h{'gene'};$temp_h{'in_dbnsfp'}=int(1); } else {$temp_h{'in_dbnsfp'}=int(0);}
   return \%temp_h
}
sub query_DBSNP {
   my %temp_h=();
   my $chr=$_[0];my $hg19pos=int($_[1]);my $alt=$_[2];my $ref=$_[3];
   my $cur_dbSNP=$DBSNP->find({'chr'=>$chr,'hg19pos'=>int($hg19pos),'alt'=>$alt,'ref'=>$ref});
   if (my $col=$cur_dbSNP->next) { %temp_h=%$col; delete $temp_h{"_id"};delete $temp_h{'chr'};delete $temp_h{'gene'};$temp_h{'in_DBSNP'}=int(1);delete $temp_h{'ref'};delete $temp_h{'alt'}} else {$temp_h{'in_DBSNP'}=int(0);}
   return \%temp_h;
}
sub CLINVAR {
   my %temp_h=();
   my $chr=$_[0];my $hg19pos=int($_[1]);my $alt=$_[2];my $ref=$_[3];
   my $cur_dbSNP=$CLINVARDB->find({'chr'=>$chr,'hg19pos'=>int($hg19pos),'alt' => $alt,'ref'=>$ref});
   if (my $col=$cur_dbSNP->next) { %temp_h=%$col; delete $temp_h{"_id"};delete $temp_h{'chr'};delete $temp_h{'gene'};delete $temp_h{'qual'};delete $temp_h{'filter'};$temp_h{'in_clinvar'}=int(1);} else {$temp_h{'in_clinvar'}=int(0);}
   return \%temp_h;
}
sub COSMIC {
   my %temp_h=();
   my $chr=$_[0];my $hg19pos=int($_[1]);my $alt=$_[2];my $ref=$_[3];
   my $cur_dbSNP=$COSMICDB->find({'chr'=>$chr,'hg19pos'=>int($hg19pos),'alt' => $alt,'ref'=>$ref});
   if (my $col=$cur_dbSNP->next) { 
     %temp_h=%$col; delete $temp_h{"_id"};delete $temp_h{'chr'};delete $temp_h{'gene'};delete $temp_h{'qual'};delete $temp_h{'filter'};$temp_h{'in_cosmic'}=int(1);
     if ($temp_h{"cosmic_info"}=~/(CNT=.*)/) { $temp_h{"cosmic_info"}=$1; }
   } else {$temp_h{'in_cosmic_exact'}=int(0);}
   return \%temp_h;
}
sub TCGA {
   my %temp_h=();
   my $chr=$_[0];my $hg19pos="$_[1]";my $alt=$_[2];my $ref=$_[3];
   $chr=~s/chr//g;
   my $cur_dbSNP=$TCGADB->find({'chr'=>"$chr",'hg19pos'=>$hg19pos,'tumor_seq_allele1' => $alt,'reference_allele'=>$ref});
   if (my $col=$cur_dbSNP->next) { %temp_h=%$col; delete $temp_h{"_id"};delete $temp_h{'hg19pos'};delete $temp_h{'chr'};delete $temp_h{'gene'};delete $temp_h{'qual'};delete $temp_h{'filter'};$temp_h{'in_tcga'}=int(1);} else {$temp_h{'in_tcga'}=int(0);}
   return \%temp_h;

}
sub SAMPLEIN {
   my %temp_h=();
   my $sample=$_[0];
   my $cur=$PATIENTS->find({'sample'=>$sample});
   if (my $col=$cur->next) { %temp_h=%$col; delete $temp_h{"_id"}; delete $temp_h{"sample"}}
   return \%temp_h;
}


sub IGENE {
   my %temp_h=();
   my $gene=$_[0];
   my $cur_gene=$GENES->find({'gene'=>$gene});
   if (my $col=$cur_gene->next) { %temp_h=%$col; delete $temp_h{"_id"};delete $temp_h{'gene'};$temp_h{'in_gene'}=int(1);} else {$temp_h{'in_gene'}=int(0);}
   return \%temp_h;
}

sub LEVEL2 {
   my %temp_h=();
   my $gene=$_[0];
   my $cur_gene=$LEVEL2G->find({'gene'=>$gene});
   if (my $col=$cur_gene->next) { %temp_h=%$col; delete $temp_h{"_id"};delete $temp_h{'gene'};$temp_h{'in_level2'}=int(1);} else {$temp_h{'in_level2'}=int(0);}
   return \%temp_h;
}

sub mongo_parse {
   my $line=$_[0];
   my @temp=split(/\s+/,$line);
   my $file_info=$_[1];
   LOOP: foreach my $set (@temp) {
     if ($set=~/=/) {
       my ($key,$val)=split("=",$set);
       $key=lc($key);
       if ($key eq "study" || $key eq "project") {next LOOP}
       $file_info->{$key}=$val;
       print OUT "$key = $val\n";
     }
   }
   print OUT "\n";
   return $file_info;
}
sub check_snpeff {
  my $var_info_ref=$_[0];
  my $key=$_[1];
  my $val=$_[2];
  my $vers=3;
  my $deg=0;
  my $store_effect="";
  my @snpeffs=split(/\,/,$val);
  my @store_temp3=();
  LOOP3: for (my $i=0;$i<=$#snpeffs;++$i) {
    if ($snpeffs[$i] =~/(.*)\((.*)\)/) {
       my $effect=$1;
       my @temp3=split(/\|/,$2);
#       if (!exists($temp3[6])) {return -1}
       if ($temp3[6] eq "protein_coding") {$vers=3}
       elsif($temp3[5] eq "protein_coding") {$vers=2}
       else {next LOOP3}
       if ($temp3[0] eq "HIGH") {
           $store_effect=$effect;
           @store_temp3=@temp3;
           $deg=4;
        } elsif ($temp3[0] eq "MODERATE" && $deg < 4) {
           $store_effect=$effect;
           @store_temp3=@temp3;
           $deg=3;
        } elsif ($temp3[0] eq "LOW" && $deg < 3) {
           $store_effect=$effect;
           $deg=2;
           @store_temp3=@temp3;
        } elsif ($temp3[0] eq "MODIFIER" && $deg < 2) {
           $store_effect=$effect;
           $deg=1;
           @store_temp3=@temp3;
        }
     }
   }
   if ($deg>2) {
      $var_info_ref->{'aberration_filter'}="PASS";
     if ($store_temp3[0] ne "") { $var_info_ref->{'functional_change'}=$store_effect}
     if ($store_temp3[2] ne "") { $var_info_ref->{'codon_change'}=$store_temp3[2];}
     if ($store_temp3[3] ne "") { $var_info_ref->{'amino_acid_change'}=$store_temp3[3];}
     if ($vers eq "2") {
       if ($store_temp3[4] ne "") { $var_info_ref->{'gene'}=$store_temp3[4]; push(@{$var_info_ref->{'genes'}},$store_temp3[4]); } 
     } elsif ($vers eq "3") {
       if ($store_temp3[5] ne "") { $var_info_ref->{'gene'}=$store_temp3[5];push(@{$var_info_ref->{'genes'}},$store_temp3[5]); }
     }
     $var_info_ref->{'snpeff'}="1";
     $var_info_ref->{'deg'}=$deg;
   } else {
      $var_info_ref->{'snpeff'}="0";
      $var_info_ref->{'aberration_filter'}="NOSNPEFF";
      return -1;
   }
   return 1;
}

