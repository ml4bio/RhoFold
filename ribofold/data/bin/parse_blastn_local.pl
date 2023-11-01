
#!/usr/bin/perl -w
use strict;



# Default options
my $E_max=10;
my $E_min=-1.0;
my $P_max=1;
my $cov_thrshd=0;
my $sc_thrshd=-10;
my $qid_thrshd=0;
my $bl=-10;            # minimum per-residue bit score with query at ends of HSP for loose end pruning
my $bs=-10;            # minimum per-residue bit score with query at ends of HSP for strict end pruning
my $bg=30;             # below this number of end gaps the loose HSP pruning score is used
my $outformat="fas";
my $append=0;
my $query_file="";
my $infile;
my $outfile;
my $v=2;

# Variable declarations
my $i;                 # residue index
my $j;                 # residue index
my $k;                 # sequence index
my $options="";
my $line;              # line read in from file
my $query_length=0;    # number of residues in query sequence
my $query_match=0;     # number of upper-case residues (=match states) in query sequence
my $capitalize=0;      # capitalize query
my $nameline;          # >template_name
my $Evalue;            # e-value of hit
my $score;             # bit score of hit
my $hit_length;        # number of residues in HSP
my $coverage;          # hit-length/$query_length
my $score_col;         # score per column
my $score_min=0;       # $score_min=-3*log($P_max)/log(2);

my $query_name;        # name of query file
my $queryseq;          # residues of query read in  with -q or -q2m option
my $qfirst;            # index of first query residue in pairwise alignment
my $qlast;             # index of last  query residue in pairwise alignment
my $tfirst;            # index of first template residue in pairwise alignment
my $tlast;             # index of last  template residue in pairwise alignment
my $tlen=0;            # length of template in pairwise alignment
my @query_res;         # query residues from current pairwise alignment
my @template_res;      # template residues from current pairwise alignment
my $query_res;         # query residues from current pairwise alignment
my $template_res;      # template residues from current pairwise alignment
my $line_number=0;
my $new_hit="";        # new sequence record; is only written if coverage threshold is exceeded
my $nhit=0;            # counts the number of sequences already in alignment
my @hitnames;          # $hitnames[$nhit] is the nameline of the ihit'th hit
my @hitseqs;           # $hitseqs[$nhit] contains the residues of the ihit'th hit
my @match;             # for -q option: $match[$i]=1 if $i'th query residue is capital letter in query, else 0
my $qid;               # $qid is sequence identity with query (for -q option: CONSIDER ONLY MATCH STATES)
my $len;               # $len is sequence number of residues of seq k aligned with a match state in query
my $b;                 # minimum per-residue bit score with query at ends of HSP
my $pfile="";          # alignment file used to calculate PSSM for -p and s/c options
my $bfile="";          # alignment file used to calculate PSSM for -b option
my $GAP=11.0/3.0;        # gap opening penalty in bits (for BLOSUM62: 11 bits/3)
my $EXTEND=1.0/3.0;      # gap extension penalty in bits (for BLOSUM62: 1 bits/3)
my @queryseq;
my $skip=0;            # skip this template sequence because it might be a synthetic fusion protein
my $best=0;            # extract only the best HSP per sequence
my $rescaled_Gonnet=0; # Gonnet matrix not yet rescaled to bits
my @qp=();             # $qb[$i][$a] is PSSM from alignment read in with -B option
my @qb=();             # $qp[$i][$a] is PSSM from alignment read in with -P option


$infile=$ARGV[0];
$query_file= $ARGV[1];
$outfile=$ARGV[2];

#Include query sequence as first sequence in alignment?
if ($query_file) {
    open(QUERYFILE,"<$query_file") or die ("ERROR: Cannot open $query_file: $!\n");
    while($line=<QUERYFILE>) # Read name line
    {
        if ($line=~/^>(.*)/)
        {
            $query_name=$1;
            last;
        }
    }
    $hitseqs[0]="";
    while($line=<QUERYFILE>)  # Read residues
    {
        if ($line=~/^>/) {last;}
        chomp($line);
        $line=~s/\s+//g;      # remove white space
        $hitseqs[0].=$line;
    }
    close(QUERYFILE);

    # Prepare name line of hit
    if ($outformat eq "psi") {
        $query_name=~/^(\S{1,20})\S*\s*(.*)/;       # delete everything after first block
        $line=sprintf("%s",$1);
        $line=~ tr/ /_/;
        $hitnames[0] = sprintf("%-31.31s ",$line);
    } else {
        $hitnames[0] = sprintf(">%s  E=0.0",$query_name);
    }
    $hitseqs[0] =~ tr/-.//d;      # delete all gaps from query
    $queryseq = $hitseqs[0];
    $hitseqs[0] =~ tr/a-z/A-Z/d;  # capitalize hitseq[0] and delete gaps
#    $hitseqs[0] =~ tr/Uu/Cc/;  # nicht mehr noetig in blast. Kann aber alignhits.pl zum abschmieren bringen.
    $nhit=1;

    # Capitalize query?
    if ($capitalize) {$queryseq =~ tr/a-z/A-Z/;}
    $query_match = ($queryseq=~tr/A-Z/A-Z/);  # count number of match states in query

    # Determine match columns as those with upper case residue in query
    @queryseq=unpack("C*",$queryseq);
    for ($j=0; $j<@queryseq; $j++) {
        if ($queryseq[$j]>=65 && $queryseq[$j]<=90) {$match[$j]=1;} else {$match[$j]=0;}
    }
}




# Scan Blast output file for query length (needed for coverage)
open(INFILE,"<$infile") or die ("Error: cannot open $infile: $!\n");
$line_number++;
while ($line=<INFILE>)
{
    if ($line=~/^Length\s*=\s*(\d+)/) {$query_length = $1; last;}
    $line_number++;
}
#print("Query length = $query_length\n");

while ($line = <INFILE>) #scan through PsiBlast-output line by line
{
    # New nameline found?
    #print "$line";
    #if ($line=~/^Length\s*=\s*(\d+)/) { print "length=$1\n\n\n\n";}

    if ($line=~s/^>//) 
    {
	#print "$line";
        $line=~s/\s+/ /g;
        chomp($line);
        $nameline=$line;
        while ($line=<INFILE>) 
	{
            if ($line=~/^Length\s*=\s*(\d+)/) {last;}
            chomp($line);
            $nameline.=$line;
        }
        $line=~/^Length\s*=\s*(\d+)/;
        $tlen=$1;
        $nameline=~s/\s+/ /g;
        $nameline=~s/\s+gi\|/   gi\|/g;
        # Is sequence a synthetic fusion protein ?
        #if ($nameline=~/(\[synthetic| synthetic|construct|cloning|vector|chimeric|fusion)/i) {$skip=1;} else {$skip=0;}

	#print "$nameline\n";
    }

    # New HSP found?
    elsif (!$skip && $line=~/^ Score =/)
    {
        if($best) {$skip=1;} # skip all following hits with same sequence?

        # First check whether E-value is small enough
        if($line =~ /^ Score =\s*(\S+)\s*bits\s*\S*\s*Expect =\s*(\S+)/) 
	{
            $score=$1;
            $Evalue=$2;

	    #print "$score, $Evalue\n";
        } 
	else 
	{
            print("\nWARNING: wrong format in blast output. Expecting Score = ... Expect = ..\n$line\n");
        }
        $Evalue=~s/^(e|E)/1$1/;   # Expect = e-123 -> 1e-123
        $Evalue=~tr/,//d;
        if ($Evalue>$E_max || $Evalue<$E_min) {$new_hit=""; next;} # reject hit

        # Record sequence identity
        # (not needed, qid calculated afterwards WITHOUT counting template residues aligned to gaps in query)
        $line=<INFILE>;
        if ($line =~ /^ Identities =\s*\S+\/(\S+)\s+\((\S+)%\)/)
	{
            $qid=$2;
	    #print "$qid\n";
            $line=<INFILE>;
        }
	else
	{
            $qid=0.0; # if match is too poor then no identities are given
        }

        # Skip another line and read following line

        $line=<INFILE>;
        $line=<INFILE>;

        # Read pairwise alignment
        $qfirst="";
        $tfirst="";
        $query_res="";
        $template_res="";
        while ($line=~/^Query\s+\d+\s+\S+\s+\d*/) # Cycle in this loop until no new "Query:" lines are found
        {	    
            if ($line!~/^Query\s+(\d+)\s+(\S+)\s+(\d*)/) 
	    {
                print("WARNING 1: wrong format of blast output in $infile, line $.\n");
                last;
            }
            if ($3 eq "") {
                <INFILE>; <INFILE>; <INFILE>; $line=<INFILE>;
                print("WARNING 2: wrong format of blast output in $infile, line $. Skipping alignment block.\n");
                next;
            }
            if ($qfirst eq "") {$qfirst=$1;}
            $query_res .= $2;
            $qlast=$3;
            <INFILE>; $line=<INFILE>;
            if ($line!~/^Sbjct\s+(\d+)\s+(\S+)\s+(\d+)/) 
	    {
                print("WARNING 3: wrong format of blast output in $infile, line $.\n");
                last;
            }
            if ($tfirst eq "") {$tfirst=$1;}
            $template_res .= $2;
            $tlast=$3;
            <INFILE>; $line=<INFILE>;
        } # end while(1)
        # Check lengths
	$query_res = uc($query_res);
        $template_res = uc($template_res);
        if (length($template_res)!=length($query_res)) {
            print("WARNING: Query and template lines do not have the same length in $infile, line $.\n");
            print("Q: $query_res\n");
            print("T: $template_res\n");
            next;
        }


	#print "$query_res\n";
	#print "$template_res\n";

        # Check whether hit has sufficient score per column
        $hit_length=($template_res=~tr/a-zA-Z/a-zA-Z/);
        if ($hit_length==0) {next;}                # Reject hit?
        $score_col=$score/$hit_length;

        @query_res =unpack("C*",$query_res);
        @template_res=unpack("C*",$template_res);

        # Prune ends of HSP which are not reliably homologous
        #if (($bs>-9 || $bl>-9) && !&PruneHSP()) {next;}   # if entire HSP is pruned away, goto next alignment

        # Check whether hit has sufficient sequence identity and coverage with query
        if (!$query_file) 
	{
            $len=0; $qid=0;
            for ($i=0; $i<scalar(@query_res); $i++) 
	    {
                if ($template_res[$i]!=45 && $query_res[$i]!=45) {  # count only non-gap template residues in match columns!
                    $len++;
                    if ($query_res[$i]==$template_res[$i]) {$qid++;}
                }
            }
            $coverage = 100*$len/$query_length;
        } 
	else 
	{
            $len=1; $qid=0; $j=$qfirst-1; # if first_res=1 then $j=0
            for ($i=0; $i<scalar(@query_res); $i++) 
	    {
                if ($query_res[$i]!=45) 
		{
                    if ($template_res[$i]!=45 && $match[$j]) {      # count only non-gap template residues in match columns!
                        $len++;
                        if ($query_res[$i]==$template_res[$i]) {$qid++;}
                    }
                    $j++;                                         # $j = next position in query
                }
            }
            $coverage = 100*$len/$query_match;
        }
        if ($len==0) {next;}                              # Reject hit?
        if (100*$qid/$len<$qid_thrshd) {next;}            # Reject hit?
        if ($coverage<$cov_thrshd) {next;}                # Reject hit?
#       print("Q: $query_res\n");
#       print("T: $template_res\n\n");

        # Check score per column
        if ($sc_thrshd>-9 || $score_min>0) {
            if (!&CheckScorePerColumn()) {next;}
        }

        if ($v>=3) {printf("nhit=%-2i  qid=%-3i  qlen=%-3i  qid=%-3i%% s/c=%-6.3f\n",$nhit,$qid,$len,100*$qid/$len,$score_col);}

        # Record residues
        $new_hit = "-"x($qfirst-1);                       # Print gaps at beginning of sequence
        if ($outformat eq "psi") {
            for ($i=0; $i<scalar(@query_res); $i++) {
                if ($query_res[$i]!=45) {                 # residues aligned to gaps are ignored
                    $new_hit .= uc(chr($template_res[$i])); # UPPER case if aligned with a query residue (match state)
                }
            }
        } else {
            for ($i=0; $i<scalar(@query_res); $i++) {
                if ($query_res[$i]!=45) {
                    $new_hit .= uc(chr($template_res[$i])); # UPPER case if aligned with a query residue (match state)
                } else {
                    if($template_res[$i]!=45) {
                        $new_hit.=lc(chr($template_res[$i])); # lower case if aligned with a gap in the query (insert state)
                    }
                }
            }
        }
        $new_hit .= "-" x ($query_length-$qlast);      # Print gaps at end of sequence
#       $new_hit =~ tr/Uu/Cc/;   # nicht mehr noetig in blast. Kann aber alignhits.pl zum abschmieren bringen.
        $hitseqs[$nhit] = $new_hit;
#       printf("%s\n",$new_hit);

        # Prepare name line of hit
        if ($outformat eq "psi") {
            $nameline=~/^(\S{1,20})\S*\s*(.*)/;           # delete everything after first block
            $line=sprintf("%s(%i-%i:%i)",$1,$tfirst,$tlast,$tlen);
            $line=~ tr/ /_/;
            $hitnames[$nhit] = sprintf("%-31.31s ",$line);
        } else {
            $nameline=~/^(\S*)\s*(.*)/;                   # delete everything after first block
            $hitnames[$nhit] = sprintf(">%s(%i-%i:%i) %s  E=%g s/c=%4.2f id=%.0f%% cov=%.0f%%",
                                       $1,$tfirst,$tlast,$tlen,$2,$Evalue,$score_col,100*$qid/$len,$coverage);
        }

        $nhit++;

	#print "$nhit\n" if($nhit%100 ==0);
    } # end elseif new HSP found
} # end while ($line)

close(INFILE);



# If output format is fasta or a2m we have to insert gaps:
if ($outformat ne "psi")
{
    my @len_ins; # $len_ins[$j] will count the maximum number of inserted residues after match state $j.
    my @inserts; # $inserts[$j] contains the insert (in small case) of sequence $k after the $j'th match state
    my $insert;
    my $ngap;

    # For each match state determine length of LONGEST insert after this match state and store in @len_ins
    for ($k=0; $k<$nhit; $k++) {
        # split into list of single match states and variable-length inserts
        # ([A-Z]|-) is the split pattern. The parenthesis indicate that split patterns are to be included as list elements
        # The '#' symbol is prepended to get rid of a perl bug in split
        $j=0;
        @inserts = split(/([A-Z]|-)/,"#".$hitseqs[$k]."#");
#       printf("%3i: %12.12s %s\n",$k,$hitnames[$k],$hitseqs[$k]);
#       printf("Sequence $k: @inserts\n");
        foreach $insert (@inserts) {
            if( !defined $len_ins[$j] || length($insert)>$len_ins[$j]) {
                $len_ins[$j]=length($insert);
            }
            $j++;
#           printf("$insert|");
        }
#       for (my $i=0; $i<@inserts; $i++) {printf("%s%-2i ",$inserts[$i],$len_ins[$i]);}
#       printf("\n");
    }

    # After each match state insert residues and fill up with gaps to $len_ins[$i] characters
    for ($k=0; $k<$nhit; $k++) {
        # split into list of single match states and variable-length inserts
        @inserts = split(/([A-Z]|-)/,"#".$hitseqs[$k]."#");
        $j=0;

        # append the missing number of gaps after each match state
        foreach $insert (@inserts) {
            if($outformat eq "fas") {
                for (my $l=length($insert); $l<$len_ins[$j]; $l++) {$insert.="-";}
            }
            else {
                for (my $l=length($insert); $l<$len_ins[$j]; $l++) {$insert.=".";}
            }
            $j++;
        }
        $hitseqs[$k] = join("",@inserts);
        $hitseqs[$k] =~ tr/\#//d; # remove the '#' symbols inserted at the beginning and end
    }
}


if ($query_file) {
    # Determine match states
    my @qa2m = unpack("C*",$hitseqs[0]); # $hitseq[0] is query sequence WITH INSERTS
    my @matchali=();
    my $L=scalar(@qa2m);
    $j=0;
    for ($i=0; $i<@match; $i++) {
        while ($j<$L && !($qa2m[$j]>=65 && $qa2m[$j]<=90)) {$matchali[$j++]=0;}  #move to column with next upper case residue
        $matchali[$j++]=$match[$i];                                           #is next query residue upper-case or not?
    }

    # Set all match states to upper case, non-match states to lower case
    my @res;
    for ($k=0; $k<$nhit; $k++) {
        @res = unpack("C*",$hitseqs[$k]);
#       printf("Q: %s\n",$hitseqs[0]);
#       printf("T: %s\n",$hitseqs[$k]);
        for ($i=0; $i<@res; $i++) {
            if ($matchali[$i]) {
                if ($res[$i]>=97 && $res[$i]<=122) {$res[$i]-=32;}  #convert to upper case
            } else {
                if ($res[$i]>=65 && $res[$i]<=90) {$res[$i]+=32;}   # convert to lower case
                elsif ($res[$i]==45) {$res[$i]=46;}     # convert '-' to '.'
            }
#           printf("%3i  Q:%s T:%s  match=%i len=%i\n",$i,chr($qa2m[$i]),chr($res[$i]),$qid[$k],$len);
        }
        $hitseqs[$k] =  pack("C*",@res);
    }
}


# Remove gaps? Captialize?
if ($outformat eq "ufas") {
    for ($k=0; $k<$nhit; $k++) {$hitseqs[$k]=~tr/a-z.-/A-Z/d;} # Transform to upper case and remove all gaps
} elsif ($outformat eq "fas") {
    for ($k=0; $k<$nhit; $k++) {$hitseqs[$k]=~tr/a-z./A-Z-/;}  # Transform to upper case
} elsif ($outformat eq "a3m") {
    for ($k=0; $k<$nhit; $k++) {$hitseqs[$k]=~tr/.//d;}        # Remove gaps aligned to inserts
}

# Write sequences into output file
open (OUTFILE, ">$outfile") or die ("cannot open $outfile:$!\n");
if ($outformat eq "psi") {
    for ($k=0; $k<$nhit; $k++) {
        $hitseqs[$k] =~ tr/./-/;
        printf(OUTFILE "%s %s\n",$hitnames[$k],$hitseqs[$k]);
    }
} 
else {
    for ($k=0; $k<$nhit; $k++) {
        printf(OUTFILE "%s\n%s\n",$hitnames[$k],$hitseqs[$k]);
    }
}
close OUTFILE;

if ($v>=1) {printf("$nhit sequences extracted from $infile and written to $outfile\n");}
exit(0);




