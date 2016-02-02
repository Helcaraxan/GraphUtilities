#!/usr/bin/perl
# sudo apt-get install perl libswitch-perl
use Switch;

my %succs=();
my %preds=();
my %index=();
my $size=0;
my $debug=0;

my %options=( "dd2" => "n", "mm" => "n", "jacobi2d" => "n T");

my $bench=shift @ARGV;
switch ($bench) {
    case "dd2"         { build_diamond2D(shift @ARGV);} # dd2 $n
    case "mm"          { build_mm(shift @ARGV);}        # mm $n
    case "jacobi2d"    { build_jacobi2d(@ARGV);}        # jacobi2d $n $T
    else { print "Usage:\n";
	   for my $b (keys %options) { print "\t./gen.pl $b ",$options{$b},"\n";}
	   die "Unknown bench $bench";
    }
}



dumpsuccs();

sub build_jacobi2d {
    my ($n,$T) = @_;

    for (my $i=0; $i<$n; $i++) {
	indx(0,$i);
    }
    for (my $t=1; $t<$T; $t++) {
	for (my $i=0; $i<$n; $i++) {
	    $preds{indx($t,$i)}{indx($t-1,$i-1)}++ if (($i-1)>=0);
	    $preds{indx($t,$i)}{indx($t-1,$i)}++;
	    $preds{indx($t,$i)}{indx($t-1,$i+1)}++ if (($i+1)<$n);
	}
    }
    preds2succs();
}

sub build_mm {
    my $n = shift @_;

    for (my $i=0; $i<$n; $i++) {
	for (my $j=0; $j<$n; $j++) {
	    for (my $k=0; $k<$n; $k++) {
		$preds{indx('C',$i,$j,$k)}{indx('C',$i,$j,$k-1)}++;
		$preds{indx('C',$i,$j,$k)}{indx('A',$i,$k)}++;
		$preds{indx('C',$i,$j,$k)}{indx('B',$k,$j)}++;
	    }
	}
    }
    preds2succs();
}

sub build_diamond2D {
    my $n = shift @_;

    for (my $i=0; $i<$n; $i++) {
	for (my $j=0; $j<$n; $j++) {
	    $preds{indx($i,$j)}{indx($i,$j-1)}++ if ($j>0);
	    $preds{indx($i,$j)}{indx($i-1,$j)}++ if ($i>0);
	}
    }
    preds2succs();
}

sub preds2succs {
    foreach my $i (keys %preds) {
	foreach my $p (keys %{$preds{$i}}) {
	    $succs{$p}{$i}++
	}
    }
}

sub dumpsuccs {
    print "$size\n";
    foreach (my $i=0;$i<$size;$i++) {
	if ($debug) {
	    print $rindx{$i}.": ", join(" ", map {$rindx{$_}} sort {$a <=> $b} keys %{$succs{$i}}), " #\n";
	} else {
	    if (scalar keys %{$succs{$i}}) {
		print "$i: ", join(" ", sort {$a <=> $b} keys %{$succs{$i}}), " #\n";
	    } else {
		print "$i: #\n";
	    }
	}
    }
}
	
sub indx {
    my $coord = join(",",@_);
    return $indx{$coord} if (defined $indx{$coord});
    $rindx{$size}="($coord)";
    $indx{$coord}=$size++;
}
    
