use warnings;
use Getopt::Long;
use List::Util qw(min max);

%len_of_chr = (
    chr1 => 249250621,
    chr2 => 243199373,
    chr3 => 198022430,
    chr4 => 191154276,
    chr5 => 180915260,
    chr6 => 171115067,
    chr7 => 159138663,
    chrX => 155270560,
    chr8 => 146364022,
    chr9 => 141213431,
    chr10 => 135534747,
    chr11 => 135006516,
    chr12 => 133851895,
    chr13 => 115169878,
    chr14 => 107349540,
    chr15 => 102531392,
    chr16 => 90354753,
    chr17 => 81195210,
    chr18 => 78077248,
    chr20 => 63025520,
    chrY => 59373566,
    chr19 => 59128983,
    chr22 => 51304566,
    chr21 => 48129895,
);

# WORKS ON 0-BASED START COORDINATES THROUGHOUT
# ASSUMES INPUT BEDGRAPH FILE SORTED AND NON-OVERLAPPING

$MIN_POS = "";
$MAX_POS = "";
$CHR = "";
$WIN_EXTEND = "";
$FEATURES_FILE = "";
GetOptions(
    "min_pos=i" => \$MIN_POS,
    "max_pos=i" => \$MAX_POS,
    "chr=s" => \$CHR,
    "features=s" => \$FEATURES_FILE,
    "win_extend=i" => \$WIN_EXTEND,
);

die if ($CHR eq "" || $FEATURES_FILE eq "" || $WIN_EXTEND eq "");
die if ($MIN_POS =~ /^\d+$/ && $MAX_POS =~ /^\d+$/ && $MIN_POS >= $MAX_POS);
if ($MIN_POS eq "") { $MIN_POS = 0; }
if ($MAX_POS eq "") { $MAX_POS = $len_of_chr{$CHR}; }
if ($FEATURES_FILE eq "stdin") {
    *FEATURES = *STDIN;
}
elsif ($FEATURES_FILE =~ /\.gz$/) {
    open FEATURES, "zcat $FEATURES_FILE | " or die $!;
}
else {
    open FEATURES, $FEATURES_FILE or die $!;
}


# Roll bedgraph file until correct chromosome.
chomp($next_line = <FEATURES>);
@F = split "\t", $next_line;
while ($F[0] ne $CHR) {
    chomp($next_line = <FEATURES>);
    @F = split "\t", $next_line;

    die if eof(FEATURES);
}

# Create first window and flush features from bedgraph file
$cur_pos = $MIN_POS;
$cur_win_start = $cur_pos - $WIN_EXTEND - 1;
$cur_win_end = $cur_pos + $WIN_EXTEND;
$cur_win_sum = $cur_win_total_olap = 0;
while ($F[2] <= $cur_win_start) {
    chomp($next_line = <FEATURES>);
    @F = split "\t", $next_line;
}

# Calculate all the overlaps for first entry
@bedgraph_entries_to_consider = ();
while ( ($olap = olap($cur_win_start, $cur_win_end, $F[1], $F[2])) > 0 ) {
    $cur_win_sum += $F[3] * $olap;
    $cur_win_total_olap += $olap;

    push @bedgraph_entries_to_consider, [$F[1], $F[2], $F[3]];  # Save the bedgraph entry

    chomp($next_line = <FEATURES>);
    @F = split "\t", $next_line;
}


# Now  @F  contains the first bedgraph entry that does not overlap with the first window anymore. 
# Print out first window, then start the loop for all subsequent windows
print_cur_window();

while ($cur_pos < $MAX_POS) {
    $prev_win_start = $cur_win_start;
    $prev_win_end = $cur_win_end;
    $cur_pos++;
    $cur_win_start = $cur_pos - $WIN_EXTEND - 1;
    $cur_win_end = $cur_pos + $WIN_EXTEND;

    # First unshift non-overlapping bedgraph entries, then push more bedgraph entries.
    # All this while keeping track of changes in sum and total_olap. 
    if (@bedgraph_entries_to_consider == 0) {
        # Do nothing.
        # No overlapping bedgraph entries with previous window
    }
    elsif ( $cur_win_start >= $bedgraph_entries_to_consider[0]->[1] ) {
        # The first bedgraph entry overlapped with previous window but does not overlap
        # with current window anymore. Unshift it. 

        $cur_win_sum -= $bedgraph_entries_to_consider[0]->[2];
        $cur_win_total_olap -= 1;
        shift @bedgraph_entries_to_consider;
    }
    else {
        # The first window overlaps with both previous window and current window.
        # Update the extend of overlap to this window. The update below takes into
        # account update in extent of overlap in both low end and high end of the window. 

        $feature_value = $bedgraph_entries_to_consider[0]->[2];
        $olap = olap($prev_win_start, $prev_win_end, $bedgraph_entries_to_consider[0]->[0], $bedgraph_entries_to_consider[0]->[1]);
        $cur_win_sum -= $feature_value * $olap;
        $cur_win_total_olap -= $olap;

        $olap = olap($cur_win_start, $cur_win_end, $bedgraph_entries_to_consider[0]->[0], $bedgraph_entries_to_consider[0]->[1]);
        $cur_win_sum += $feature_value * $olap;
        $cur_win_total_olap += $olap;
    }

    # Bring a new bedgraph entries in if possible. Can get there only if chromosomes match. 
    # So if more data to come in bedgraph file but are in next chromosome, then will be stalled here. 
    if ( $CHR eq $F[0] && $cur_win_end > $F[1] ) {
        $cur_win_sum += $F[3];
        $cur_win_total_olap += 1;

        push @bedgraph_entries_to_consider, [$F[1], $F[2], $F[3]];  # Save the bedgraph entry

        # Start the next unoverlapped bedgraph entry, if more stuff to come from file
        if (eof(FEATURES)) {
            # If reached end of bedgraph file, just reset bedgraph entry data so
            # that nothing more is read. 

            $F[0] = "";
        }
        else {
            chomp($next_line = <FEATURES>);
            @F = split "\t", $next_line;
        }
    }
    elsif (@bedgraph_entries_to_consider > 1 && $bedgraph_entries_to_consider[$#bedgraph_entries_to_consider]->[1] > $prev_win_end) {
        # Else if the rightmost bedgraph entry not different to the leftmost, then
        # update the current overlap to the rightmost entry if possible

        $rightmost_idx = $#bedgraph_entries_to_consider;
        $cur_win_sum += $bedgraph_entries_to_consider[$rightmost_idx]->[2];
        $cur_win_total_olap += 1;
    }

    print_cur_window();
}



sub olap {
    if ($_[0] < $_[3] && $_[2] < $_[1]) {
        return( min($_[1], $_[3]) - max($_[0], $_[2]) );
    }
    else {
        return -1;
    }
}

sub print_cur_window {
    if ($cur_win_total_olap == 0) {
        print "$cur_pos\tNA\n";
    }
    else {
        $rounded = sprintf "%.2f", $cur_win_sum/$cur_win_total_olap;
        print("$cur_pos\t", $rounded, "\n");
    }
}
