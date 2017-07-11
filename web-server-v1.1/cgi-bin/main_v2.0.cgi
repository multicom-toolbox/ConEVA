#!/usr/bin/perl
# 12/25/2015, Badri Adhikari

use strict;
use warnings;
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser); 
use File::Basename;
use lib '/var/www/cgi-bin/coneva/Roman-1.24/lib';
use Roman;
use Time::Piece;

$CGI::POST_MAX = 1024 * 1000000; # 100MB
$CGI::DISABLE_UPLOADS = 0;

# Input Parameters
my $query = new CGI;
my $param_cont_type  = $query -> param("contacttype");
my $param_rr_raw     = $query -> param("rr_raw");        # this can be a contacts zip file as well
my $param_rr_file    = $query -> param("rrfile");
my $param_na_rr_fl   = $query -> param("native_rr_flag");
my $param_pdb_raw    = $query -> param("pdb_raw");
my $param_pdb_file   = $query -> param("pdbfile");
my $param_sec        = $query -> param("protein_sec");
my $param_atom_type  = $query -> param("atomtype");      # ca/cb/any/heavyatoms
my $d_threshold      = $query -> param("distthreshold");
my $short_range_min  = $query -> param("srmin");         # 6;
my $short_range_max  = $query -> param("srmax");         # 11;
my $medum_range_min  = $query -> param("mrmin");         # 12;
my $medium_range_max = $query -> param("mrmax");         # 23;
my $long_range_min   = $query -> param("lrmin");         # 24;
my $long_range_max   = $query -> param("lrmax");         # 100000;
my $param_max_con    = $query -> param("maxcon");        # L/2L/ALL - Number of contacts to load from each source
my $JN               = $query -> param("neighborsize");  # Neighborhood size for computing Jaccard similarity matrix
$JN                  = int($JN);
my $MAX_DIAGRAMS     = 3;

my $remote_addr    = $ENV{'REMOTE_ADDR'};
my $now = (localtime); $now =~ s/\s+/_/g; $now =~ s/:/_/g;
error_exit("ERROR! Invalid IP address!") if length($remote_addr) < 5;

# http://proteopedia.org/wiki/index.php/Standard_Residues and http://prody.csb.pitt.edu/manual/reference/atomic/flags.html
our %AA3TO1 = qw(ALA A ASN N CYS C GLN Q HIS H LEU L MET M PRO P THR T TYR Y ARG R ASP D GLU E GLY G ILE I LYS K PHE F SER S TRP W VAL V UNK -);
our %AA1TO3 = reverse %AA3TO1;

my $APP_ROOT = "/var/www/cgi-bin/coneva";
my $program_dssp = "$APP_ROOT/dssp-2.0.4-linux-amd64";

my $html_display_width = 1000;

# Global Variables
my $MAX_SEQ_LIMIT = 1000;
my $MIN_SEQ_LIMIT = 5;
my $MAX_SOURCES   = 10;
my $MAX_JOBS_PER_CLIENT = 25;

# Only predicted top "$MAX_CONTACTS times L" contacts are read, for fast processing
# Later, I find that it is better to remove this features as it (2) adds complexity to understand, and
# (2) creates a problem when computing AUC and MCC.
my $MAX_CONTACTS  = 100; # Load 100L (all) contacts by default
#$MAX_CONTACTS = 1 if $param_max_con eq "L";
#$MAX_CONTACTS = 2 if $param_max_con eq "2L";

my $rr_sequence   = undef;
my $fasta         = undef;
my $native_pdb    = undef;
my $native_rr     = undef;
my $input_rr      = undef;
my $input_sec     = undef;
my $flag_1_src    = 0;
my $safe_filename_chars = "a-zA-Z0-9_.-";
my %top_count      = ();
my %top_count_order= ();
my %top_count_desc = ();
my %sources        = ();
my $cont_type_desc = ();
my %coord_number   = ();
my %precision      = ();
my %coverage       = ();
my %mcc            = ();
my %xd             = ();
my %spread         = ();
my %fp             = ();
my %min_true_dist  = ();
my %min_true_d_atoms = ();
my %true_contacts  = ();
my %native_res_list = ();
my %selected_contacts = ();
my $total_input_rr = undef;
my $ignored_rr_cnt = undef;
my $sectionID      = 1;

# Assess All contacts or Long-Range only
if ($param_cont_type eq "all"){
	$cont_type_desc = "ALL";
}
if ($param_cont_type eq "long"){
	$cont_type_desc = "Long-Range";
}
if ($param_cont_type eq "short"){
	$cont_type_desc = "Short-Range";
}
if ($param_cont_type eq "medium"){
	$cont_type_desc = "Medium-Range";
}

# The overall min_seq_sep for all analysis
my $min_seq_sep = $short_range_min;
$min_seq_sep = $medum_range_min if $param_cont_type eq "medium";
$min_seq_sep = $long_range_min if $param_cont_type eq "long";

# Start showing something to the browser
print "Content-type: text/html\n\n";
open HEADER, "/var/www/html/coneva/header.php"; 
while(<HEADER>){
	print $_;
}
close HEADER;

print "\n<div id=\"progressbardiv\" width=\"100px\" align=\"center\" class=\"demo-wrapper html5-progress-bar\">";
print "\n<div class=\"progress-bar-wrapper\">";
print "\n<progress id=\"progressbar\" value=\"0\" max=\"100\"></progress>";
print "\n</div>";
print "\n</br>";
print "\n</div>";
print "\n<script type=\"text/javascript\">updateProgress(\"1\");</script>";

my $file_history = "/var/www/html/coneva/jobs/history.log";
system("rm -f /var/www/html/coneva/jobs/clientmap");
my $file_clientmap = "/var/www/html/coneva/jobs/map";
my %client_map;
open CLIENTS, $file_clientmap; 
while(<CLIENTS>){
	my @C = split /\s+/, $_;
	$client_map{$C[0]} = $C[1];
}
close CLIENTS;

if(not defined $client_map{$remote_addr}){
	my $datestring = localtime->strftime('%Y-%m-%d-%H-%M-%S');
	print2file($file_clientmap, "$remote_addr $datestring");
	$client_map{$remote_addr} = $datestring;
}
$remote_addr = $client_map{$remote_addr};

# Create 1 directory for each client
my $dir_client = "/var/www/html/coneva/jobs/$remote_addr";
system("mkdir -p $dir_client");
# Recycle job folders if the client submits a lot of jobs
my $job_id = "job_id_1";
if (not -f "$dir_client/most_recent_job_id.txt"){
	print2file($file_history, "");
	print2file($file_history, "Creating first job for client $remote_addr");
	system("echo 1 > $dir_client/most_recent_job_id.txt");
}
else{
	my $current_job_id = `cat $dir_client/most_recent_job_id.txt`;
	chomp $current_job_id;
	if($current_job_id >= $MAX_JOBS_PER_CLIENT){
		$current_job_id = 0;
		print2file($file_history, "Job count limit reached for client $remote_addr");
	}
	$job_id = $current_job_id + 1;
	system("echo $job_id > $dir_client/most_recent_job_id.txt");
	print2file($file_history, "");
	print2file($file_history, "Creating job_id $job_id for client $remote_addr");
	$job_id = "job_id_$job_id";
}
# Make folders and log files for this job
system("mkdir -p $dir_client/$job_id");
system("rm -rf $dir_client/$job_id/*");
my $dir_job = "$dir_client/$job_id";
chdir $dir_job or confess "Could not change to $dir_job";
my $file_log = "$dir_job/log.txt";
mkdir "$dir_job/chimera-scripts" or confess $!;
mkdir "$dir_job/not-in-native" or confess $!;
mkdir "$dir_job/roc" or confess $!;
mkdir "$dir_job/auc_pr" or confess $!;
mkdir "$dir_job/chord" or confess $!;
mkdir "$dir_job/cmap" or confess $!;
mkdir "$dir_job/jaccard" or confess $!;
my $web_root = "http://iris.rnet.missouri.edu/coneva/jobs/$remote_addr/$job_id";

# Create input files
chk_errors_and_log_params();
create_files();
chk_errors_input_files();

# Precompute this once for all, to save time
if ($native_pdb){
	my $min_seq_sep = $short_range_min;
	$min_seq_sep = $medum_range_min if $param_cont_type eq "medium";
	$min_seq_sep = $long_range_min if $param_cont_type eq "long";
	%min_true_dist = all_pairs_min_dist($native_pdb, $min_seq_sep, $d_threshold, 0);
	%min_true_d_atoms = all_pairs_min_dist($native_pdb, $min_seq_sep, $d_threshold, 1);
	pdb2rr($native_pdb, "native.rr", $min_seq_sep, $d_threshold, $param_atom_type);
	$native_rr = "native.rr";
	%true_contacts = rrfile_to_r1r2hash($native_rr, $min_seq_sep, 1000000);
	%native_res_list = res_num_res_name($native_pdb);
}

# Obtain EVAcon format native RR file
if ($native_pdb){
	print2file("native-evacon.rr", "FRMAT RR");
	print2file("native-evacon.rr", "TARGET xxxxx");
	print2file("native-evacon.rr", "AUTHOR Author's name");
	print2file("native-evacon.rr", "REMARK Predictor remarks");
	print2file("native-evacon.rr", "METHOD Description of methods used");
	my $natseq = seq_chain_with_gaps($native_pdb);
	for(my $i=0; $i <= length($natseq); $i = $i + 50){
		print2file("native-evacon.rr", "SEQRES ".(substr $natseq, $i, 50));
	}
	print2file("native-evacon.rr", "MODEL 1");
	foreach (sort keys %true_contacts){
		my @C = split /\s+/, $_;
		print2file("native-evacon.rr", "CONTC ".$C[0]." ".$C[1]." ".$AA3TO1{$native_res_list{$C[0]}}." ".$AA3TO1{$native_res_list{$C[1]}}." 0 ".$d_threshold." ".$true_contacts{$_});
	}
	print2file("native-evacon.rr", "END");
}

# Kinds of inputs:
# Case 1: No native: 1D and 2D visualizations only
# Case 2: Native only
# Case 3: Native with one set of contacts
# Case 4: Native with more than one set of contacts, i.e. multiple sources
my $input_case    = undef;
if (not -f $native_pdb){
	$input_case = 1;
}
elsif ($param_na_rr_fl eq "on"){
	$input_case = 2;
}
elsif (-f $native_pdb and $flag_1_src){
	$input_case = 3;
}
elsif (-f $native_pdb and !$flag_1_src){
	$input_case = 4;
}
else{
	confess "Case not defined" if not defined $input_case;
}

if ($flag_1_src or $param_na_rr_fl eq "on"){
	$MAX_DIAGRAMS = 6;
}
if ($input_case >= 2){
	chk_rr_not_in_native();
}

load_top_count($input_rr);
load_sources($input_rr);

print_overview();
set_progress(1);
print_contact_counts();
set_progress(5);
if (!($flag_1_src or $param_na_rr_fl eq "on")){
	if ((scalar keys %sources) > 2){
		print_jaccard_similarity();
	}
}
set_progress(15);
print_chord_diagram();
set_progress(25);

if($input_case == 1){
	print_contact_maps();
	set_progress(50);
	print_coord_num_rows();
	set_progress(100);
}
elsif($input_case == 2){
	print_contact_maps();
	set_progress(30);
	print_coverage();
	set_progress(35);
	print_xd();
	set_progress(50);
	print_spread();
	set_progress(80);
	print_coord_num_rows();
	set_progress(100);
}
elsif($input_case == 3){
	print_contact_maps();
	set_progress(30);
	# Printing the native cmap separately, because it makes no sense to compare "predicted top x" and "native random x"
#	print_native_contact_maps();
	set_progress(35);
	# Print precision related things only for input.rr (because native has all 100%)
	print_precision();
	set_progress(40);
	print_FP_assessment();
	set_progress(45);
	# Print other metrics for both native.rr and input.rr
	mkdir "input" or confess $!;
	system_cmd("cp native.rr ./input/native");
	system_cmd("cp input.rr ./input/input");
	system_cmd("mv input.rr input_original.rr");
	make_rr_from_input_folder("input");
	$input_rr = "$dir_job/input.rr";
	load_sources($input_rr);
	load_top_count($input_rr);
	print_coverage();
	set_progress(50);
	print_xd();
	set_progress(55);
	print_spread();
	set_progress(60);
	# Again, use the input contacts only
	$input_rr = "input_original.rr";
	load_sources($input_rr);
	load_top_count($input_rr);
	$input_rr = "$dir_job/input.rr";
	load_sources($input_rr);
	load_top_count($input_rr);
	print_roc();
	set_progress(70);
	print_auc_pr();
	set_progress(80);
	print_MCC();
	set_progress(90);
	print_coord_num_rows();
	set_progress(100);
}
elsif($input_case == 4){
	print_contact_maps();
	set_progress(30);
	# Printing the native cmap separately, because it makes no sense to compare "predicted top x" and "native random x"
#	print_native_contact_maps();
	set_progress(35);
	print_precision();
	set_progress(40);
	print_FP_assessment();
	set_progress(45);
	print_coverage();
	set_progress(50);
	print_xd();
	set_progress(55);
	print_spread();
	set_progress(60);
	print_roc();
	set_progress(70);
	print_auc_pr();
	set_progress(80);
	print_MCC();
	set_progress(90);
	print_coord_num_rows();
	set_progress(100);
}

print_footer();
exit;

sub print_overview{
	print "\n<table class=\"results\">";
	print "\n<th colspan=2 class=\"resultsheadingrow\"><b>Section ".Roman($sectionID++)."</b> Summary <font class=\"showdescription\"><a href=\"$web_root\" target=\"_blank\">job-dir</a></font>
	</th>";
	if (defined $param_rr_raw and $param_rr_raw =~ /^http/){
		print "\n<tr><td>RR Input</td><td class=\"resultdata\">$param_rr_raw</td></tr>";
	}
	if (defined $param_rr_file){
		print "\n<tr><td>RR Input</td><td class=\"resultdata\">$param_rr_file</td></tr>";
	}
	if (defined $param_pdb_raw and $param_pdb_raw =~ /^http/){
		print "\n<tr><td>PDB Input</td><td class=\"resultdata\">$param_pdb_raw  <a target=\"_blank\" href=\"$web_root/native-evacon.rr\">computed RR file</a> </td></tr>";
	}
	if (defined $param_pdb_file){
		print "\n<tr><td>PDB Input</td><td class=\"resultdata\">$param_pdb_file  <a target=\"_blank\" href=\"$web_root/native-evacon.rr\">computed RR file</a> </td></tr>";
	}
	print "\n<tr><td>Sequence Length</td>";
	print "\n<td class=\"resultdata\">";
	print "\n RR - ".(length $rr_sequence);
	if (-f $native_pdb){
		my $seq = seq_chain_with_gaps($native_pdb);
		print "\n, Native - ".(length $seq);
		my $warning = "";
		if ((length $seq) != (length $rr_sequence)){
			$warning = "<font color=\"red\">!</font>";
		}
		print "\n <a target=\"_blank\" href=\"$web_root/sequence_compare.html\">sequence comparison</a>$warning";
	}
	print "\n</td></tr>";
	#my $maxconinfo = $param_max_con;
	#$maxconinfo = "Top-$maxconinfo" if $param_max_con ne "ALL";
	#print "\n<tr><td>Contact Pre-selection</td><td class=\"resultdata\">$maxconinfo contacts are selected from each RR file";
	#print "\n</td></tr>";
	my $atom = $param_atom_type;
	$atom = "C&alpha;" if $param_atom_type eq "ca";
	$atom = "C&beta;" if $param_atom_type eq "cb";
	print "\n<tr><td>Contact Type</td><td class=\"resultdata\">$cont_type_desc contacts with $atom atoms within $d_threshold &#197 with minimum $min_seq_sep residues sequence separation";
	print "\n</td></tr>";
	print "\n</table></td></tr></table></br>";
}

sub print_contact_counts{
	print "\n<table class=\"results\">";
	print "\n<th colspan=100 class=\"resultsheadingrow\"><b>Section ".Roman($sectionID++)."</b> Contact Count <font class=\"showdescription\"><a onclick=\"toggle_it('contactcountdesc')\">see description</a></font></th>";
	# Hidden description
	print "\n<tr style=\"display:none\" id=\"contactcountdesc\"><td colspan=100 class=\"resultdata\"><ul>";
	print "\n<li> <b>RR-File</b> column lists all the sources(sets) of input contacts. If only one set of contacts is supplied it is referred <b>\"all\"</b>.</li>";
	print "\n<li> <b>Total</b> column shows the maxium number of contacts(that agree the selected contact definition) received for that source. If <b>\"Max. No. of Top Ranked Contacts to load from each RR file\"</b> parameter is selected as <b>all</b>, all the contacts from each source are processed. However, predicted contacts that have very small sequence separation, i.e. sequence separation smaller than the definition of short-range contacts are ignored.</li>";
	print "\n<li> <b>Not in Native</b> column shows the count of the ones not found in native. In the native PDB (if supplied) if any of the input contacts are not defined, they are removed for calculations and visualizations. For instance, some residues may not be defined in the PDB file [column 2].</li>";
	print "\n<li> <b>Remaining Total</b> column shows actual number of contacts used for analysis,  evaluation and visualization.</li>";
	print "\n<li> <b>Short-range</b>, <b>Medium-Range</b>, and <b>Long-Range</b> columns show the proportion of short-, medium-, and long-range contacts in the <b>Remaining Total</b> contacts. The '.com' UCSF Chimera scripts can be downloaded by clicking the numbers for visualization.</li>";
	print "\n<li>The last six columns show the number of contacts in the <b>top-5</b>, <b>top-L/10</b>, <b>top-L/5</b>, <b>top-L/2</b>, <b>top-L</b>, and <b>top-2L</b> subsets. When native contacts are being analyzed, these top subsets are random subsets.</li>";
	print "\n</ul></td></tr>";
	# Results
	print "\n<tr>";
	my $note = "L (= ".length(seq_rr($input_rr)).") from Input RR (ignoring gaps)";
	$note = "L (= ".length(seq_chain($native_pdb)).") from Native PDB (ignoring gaps)" if -f $native_pdb;
	print "\n<td class=\"contactcountheader\" rowspan=2>&nbsp;RR-File&nbsp;</td>";
	print "\n<td class=\"contactcountheader\" colspan=6>$cont_type_desc Contacts</td>";
	print "\n<td class=\"contactcountheader\" colspan=6>$note</td>";
	print "\n</tr>";
	print "\n<tr>";
	print "\n<td class=\"contactcountheader\">&nbsp;Total&nbsp;</td>";
	print "\n<td class=\"contactcountheader\">&nbsp;Not in</br>Native&nbsp;</td>";
	print "\n<td class=\"contactcountheader\">&nbsp;Remaining Total&nbsp;</td>";
	print "\n<td class=\"contactcountheader\">&nbsp;Short-Range&nbsp;</td>";
	print "\n<td class=\"contactcountheader\">&nbsp;Medium-Range&nbsp;</td>";
	print "\n<td class=\"contactcountheader\">&nbsp;Long-Range&nbsp;</td>";
	foreach my $top (sort {$top_count{$a} <=> $top_count{$b}} keys %top_count){
		$top =~ s/L/L\// if $top =~ /L[0-9]/;
		my $prefix = "";
		$prefix = "Random" if $param_na_rr_fl;
		print "\n<td class=\"contactcountheader\">&nbsp;${prefix} $top&nbsp;</td>";
	}
	print "\n</tr>";
	foreach my $src (sort keys %sources){
		print "\n<tr>";
		# File Name
		print "\n<td class=\"contactcountvalues\">$src&nbsp;</td>";
		my %rr = rr_rows_ordered_in_hash($input_rr, 1000000, $src, "all");
		my $ignored_rr_cnt = 0;
		if (-f "./not-in-native/not_in_native_$src.rr"){
			my %rr_ignored = rr_rows_ordered_in_hash("./not-in-native/not_in_native_$src.rr", 1000000, $src, "all") ;
			$ignored_rr_cnt = scalar keys %rr_ignored;
		}
		$total_input_rr = $ignored_rr_cnt + (scalar keys %rr);
		system_cmd("rm -f temp.rr");
		print2file("temp.rr", seq_rr($input_rr));
		foreach (sort {$a <=> $b} keys %rr){
			print2file("temp.rr", $rr{$_});
		}
		# Total contacts
		print "\n<td class=\"contactcountvalues\">$total_input_rr</td>";
		# Contacts not in PDB (ignored)
		if ($ignored_rr_cnt){
			print "\n<td class=\"contactcountvalues\"><a style=\"text-decoration: none; color:#000099;\" target=\"_blank\" href=\"$web_root/not-in-native/not_in_native_$src.rr\">".$ignored_rr_cnt."</a>&nbsp;</td>";
		}
		else{
			print "\n<td class=\"contactcountvalues\">"."-"."&nbsp;</td>";
		}
		my %sr  = rr_rows_ordered_in_hash("temp.rr", 1000000, "all", "short");
		my %mr  = rr_rows_ordered_in_hash("temp.rr", 1000000, "all", "medium");
		my %lr  = rr_rows_ordered_in_hash("temp.rr", 1000000, "all", "long");
		my %all = rr_rows_ordered_in_hash("temp.rr", 1000000, "all", "all");
		my $srcnt;
		my $mrcnt;
		my $lrcnt;
		my $alcnt;
		if (-f $native_pdb){
			my %temp = rrhash2dist(\%sr, \%min_true_dist);
			foreach(keys %sr){
				my @C = split /\s+/, $sr{$_};
				delete $sr{$_} if not defined $temp{$sr{$_}};
			}
			rrhash2chimerascript(\%sr, "./chimera-scripts/${src}_sr.com");
			%temp = rrhash2dist(\%mr, \%min_true_dist);
			foreach(keys %mr){
				my @C = split /\s+/, $mr{$_};
				delete $mr{$_} if not defined $temp{$mr{$_}};
			}
			rrhash2chimerascript(\%mr, "./chimera-scripts/${src}_mr.com");
			%temp = rrhash2dist(\%lr, \%min_true_dist);
			foreach(keys %lr){
				my @C = split /\s+/, $lr{$_};
				delete $lr{$_} if not defined $temp{$lr{$_}};
			}
			rrhash2chimerascript(\%lr, "./chimera-scripts/${src}_lr.com");
			%temp = rrhash2dist(\%all, \%min_true_dist);
			foreach(keys %all){
				my @C = split /\s+/, $all{$_};
				delete $all{$_} if not defined $temp{$all{$_}};
			}
			rrhash2chimerascript(\%all, "./chimera-scripts/${src}_all.com");
		}
		$srcnt = scalar keys %sr;
		$mrcnt = scalar keys %mr;
		$lrcnt = scalar keys %lr;
		$alcnt = scalar keys %all;
		if(-f $native_pdb and $alcnt){
			print "\n<td class=\"contactcountvalues\"><a href=\"$web_root/chimera-scripts/${src}_all.com\">$alcnt</a></td>";
		}
		elsif($alcnt){
			print "\n<td class=\"contactcountvalues\">$alcnt</td>";
		}
		else{
			print "\n<td class=\"contactcountvalues\">-</td>";
		}
		if(-f $native_pdb and $srcnt){
			print "\n<td class=\"contactcountvalues\"><a href=\"$web_root/chimera-scripts/${src}_sr.com\">$srcnt<font style=\"font-family:monospace\">(".(sprintf "%.0f", 100*$srcnt/$alcnt)."%)</font></a></td>";
		}
		elsif($srcnt){
			print "\n<td class=\"contactcountvalues\">$srcnt<font style=\"font-family:monospace\">(".(sprintf "%.0f", 100*$srcnt/$alcnt)."%)</td>";
		}
		else{
			print "\n<td class=\"contactcountvalues\">-</td>";
		}
		if(-f $native_pdb and $mrcnt){
			print "\n<td class=\"contactcountvalues\"><a href=\"$web_root/chimera-scripts/${src}_mr.com\">$mrcnt<font style=\"font-family:monospace\">(".(sprintf "%.0f", 100*$mrcnt/$alcnt)."%)</font></a></td>";
		}
		elsif($mrcnt){
			print "\n<td class=\"contactcountvalues\">$mrcnt<font style=\"font-family:monospace\">(".(sprintf "%.0f", 100*$mrcnt/$alcnt)."%)</font></td>";
		}
		else{
			print "\n<td class=\"contactcountvalues\">-</td>";
		}
		if(-f $native_pdb and $lrcnt){
			print "\n<td class=\"contactcountvalues\"><a href=\"$web_root/chimera-scripts/${src}_lr.com\">$lrcnt<font style=\"font-family:monospace\">(".(sprintf "%.0f", 100*$lrcnt/$alcnt)."%)</font></a></td>";
		}
		elsif($lrcnt){
			print "\n<td class=\"contactcountvalues\">$lrcnt<font style=\"font-family:monospace\">(".(sprintf "%.0f", 100*$lrcnt/$alcnt)."%)</td>";
		}
		else{
			print "\n<td class=\"contactcountvalues\">-</td>";
		}
		foreach my $top (sort {$top_count{$a} <=> $top_count{$b}} keys %top_count){
			my %temp = rr_rows_ordered_in_hash("temp.rr", $top_count{$top}, "all", $param_cont_type);
			if(-f $native_pdb){
				rrhash2chimerascript(\%temp, "./chimera-scripts/${src}_$top.com");
				print "\n<td class=\"contactcountvalues\"><a href=\"$web_root/chimera-scripts/${src}_$top.com\" target=\"_blank\">".(scalar keys %temp)."</a></td>";
			}
			else{
				print "\n<td class=\"contactcountvalues\">".(scalar keys %temp)."</td>";
			}
		}
		print "\n</tr>";
	}
	if(-f $native_pdb){
		print "\n<tr><td colspan=13 align=left style=\"color:blue;\"></br><u><a onclick=\"toggle_it('chimeravisualize')\">How to visualize these contacts in UCSF Chimera?</a></u></td></tr>";
	}
	print "\n<tr id=\"chimeravisualize\" style=\"display:none\"><td colspan=13><font style=\"font-family:monospace\">";
	print "\n&nbsp; 1. Download and install UCSF Chimera 1.10 or later (if you do not have it already).</br>";
	print "\n&nbsp; 2. Right click a number above and select <b>'Save link as...'</b> to download a UCSF Chimera script. The scripts have '.com' extension.</br>";
	print "\n&nbsp; 3. In UCSF Chimera, use <b>File->Open</b> to open the downloaded '.com' commands script file.</br>";
	print "\n&nbsp; Chimera will automatically retrieve and load the native structure and show predicted contacts within the native residues along with their true distances.</br>";
	print "\n</font></td></tr>";
	print "\n</table></br>";
}

sub print_jaccard_similarity{
	append2log("Calculate and plot Jaccard similarity matrices..");
	chdir "$dir_job/jaccard" or confess $!;
	my $i = 1;
	foreach my $top (sort {$top_count_order{$a} <=> $top_count_order{$b}} keys %top_count_order){
		append2log("Calculate for Top $top");
		my $prefix = "";
		$prefix = "Random-" if $param_na_rr_fl eq "on";
		my $id = "${prefix}Top-$top";
		$id =~ s/L/L\// if $id =~ /L[0-9]/;
		my $out = "mat_top$top.txt";
		foreach my $src (sort keys %sources){
			my $name = $src;
			$name = substr($src, 0, 4)."..".substr($src, length($src)-4, 4) if length($src) > 8;
			print2file($out, "\t".$name, "");
		}
		print2file($out, "");
		foreach my $src1 (sort keys %sources){
			my %rr_src1 = %{$selected_contacts{$src1." ".$top}};
			confess "\nERROR! Empty RR file! Not expecting!" if not scalar keys %rr_src1;
			my $name = $src1;
			$name = substr($src1, 0, 4)."..".substr($src1, length($src1)-4, 4) if length($src1) > 8;
			print2file($out, $name."\t", "");
			foreach my $src2 (sort keys %sources){
				my %rr_src2 = %{$selected_contacts{$src2." ".$top}};
				confess "\nERROR! Empty RR file! Not expecting!" if not scalar keys %rr_src2;
				print2file($out, calc_jaccard_similarity($src1, $src2, \%rr_src1, \%rr_src2)."\t", "");
			}
			print2file($out, "");
		}
		append2log("Plot for Top $top..");
		print2file("jc$top.R", "if (!require(\"gplots\")) {");
		print2file("jc$top.R", "  install.packages(\"gplots\", dependencies = TRUE)");
		print2file("jc$top.R", "  library(gplots)");
		print2file("jc$top.R", "}");
		print2file("jc$top.R", "if (!require(\"RColorBrewer\")) {");
		print2file("jc$top.R", "  install.packages(\"RColorBrewer\", dependencies = TRUE)");
		print2file("jc$top.R", "  library(RColorBrewer)");
		print2file("jc$top.R", "}");
		print2file("jc$top.R", "png(\"$top.png\", width=2000, height=2000, res=300)");
		print2file("jc$top.R", "table_jaccard = read.table(file=\"mat_top${top}.txt\")");
		print2file("jc$top.R", "jacc_matrix = data.matrix(table_jaccard)");
		print2file("jc$top.R", "heatmap.2(jacc_matrix, tracecol = NA, symm = T, dendrogram = \"col\", Colv = T, scale = 'none', col=colorRampPalette(brewer.pal(8,\"YlOrRd\"))(length(jacc_matrix)), key.ylab = NA, key.xlab = \"Jaccard Similarity\", lmat = rbind(c(3,4),c(1,2)), lwid = c(4,1), lhei = c(1,4), density.info=\"none\", margins = c(12, 12), srtCol=90, main=\"$id\")");
		system_cmd("Rscript --no-save --no-restore  --verbose jc$top.R &> Rscript.log");
		#system_cmd("Rscript --no-save --no-restore --verbose jc$top.R 2>&1"); # for debugging
		$i++;
		last if $i == $MAX_DIAGRAMS + 1;
	}
	append2log("Printing Jaccard Similarity");
	print "\n<table class=\"results\">";
	print "\n<th colspan=13 class=\"resultsheadingrow\"><b>Section ".Roman($sectionID++)."</b> Jaccard Similariy Matrices with Dendrograms [N=$JN] <font class=\"showdescription\"><a onclick=\"toggle_it('jaccdesc')\">see description</a></font></th>";
	# Hidden description
	print "\n<tr style=\"display:none\" id=\"jaccdesc\"><td colspan=100 class=\"resultdata\"><ul>";
	print "\n<li>The Jaccard coefficient measures similarity between two sets of contacts. It is defined as the size of the intersection divided by the size of the union of the two sets.</li>";
	print "\n<li>Jaccard Similarity = |A &cap; B| / |A &cup; B|.</li>";
	print "\n<li>When N = 0, contacts are not relaxed.</li>";
	print "\n<li>When relaxing the similarity computation using neiborhood contacts, all contacts with ±1 residue deviation are considered as 1 contact. For instance, the contacts 12-30, 12-31, 11-31, 11-30, 12-31, etc. will be considered as 1 contact while computing the intersections and unions.</li>";
	print "\n</ul></td></tr>";
	print "\n<tr>";
	# Results
	my $i = 1;
	foreach my $top (sort {$top_count{$a} <=> $top_count{$b}} keys %top_count){
		my $prefix = "";
		$prefix = "Random-" if $param_na_rr_fl eq "on";
		my $id = "${prefix}Top-$top";
		$id =~ s/L/L\// if $id =~ /L[0-9]/;
		print "\n<td align=\"center\">";
		print "\n<a href=\"$web_root/jaccard/$top.png\" target=\"_blank\"><img src=\"$web_root/jaccard/$top.png\" width=\"".(2 * $html_display_width / (scalar keys %top_count) - 5)."\" align=\"middle\"/></a>";
		print "\n</br><a href=\"$web_root/jaccard/mat_top$top.txt\" target=\"_blank\">see matrix data</a>";
		if ($i == int((scalar keys %top_count)/2)){
			print "\n</td></tr><tr>";
		}
		else{
			print "\n</td>";
		}
		$i++;
		last if $i == $MAX_DIAGRAMS + 1;
	}
	print "\n</tr>";
	print "\n</table>";
	print "\n</br>";
	chdir $dir_job or confess $!;
}

sub print_chord_diagram{
	append2log("Print Chord diagram..");
	chdir "$dir_job/chord" or confess $!;
	my $i = 1;
	foreach my $top (sort {$top_count_order{$a} <=> $top_count_order{$b}} keys %top_count_order){
		append2log("Calculate for Top $top");
		my $prefix = "";
		$prefix = "Random-" if $param_na_rr_fl eq "on";
		my $id = "${prefix}Top-$top";
		$id =~ s/L/L\// if $id =~ /L[0-9]/;
		my $out = "input_top$top.rr";
		print2file($out, $rr_sequence);
		foreach my $src (sort keys %sources){
			my %rows = rr_rows_ordered_in_hash($input_rr, $top_count{$top}, $src, $param_cont_type);
			foreach (keys %rows){
				print2file($out, $rows{$_});
			}
		}
		append2log("Generate R script for for Top $top..");
		system_cmd("python $APP_ROOT/chord_diagram.py input_top$top.rr top$top.R top$top.png $id");
		append2log("Plot for Top $top..");
		#sleep 1;
		system("Rscript top$top.R");
		$i++;
		last if $i == $MAX_DIAGRAMS + 1;
	}
	print "\n<table class=\"results\">";
	print "\n<th colspan=13 class=\"resultsheadingrow\"><b>Section ".Roman($sectionID++)."</b> Chord Diagrams <font class=\"showdescription\"><a onclick=\"toggle_it('chorddesc')\">see description</a></font></th>";
	# Hidden description
	print "\n<tr style=\"display:none\" id=\"chorddesc\"><td colspan=100 class=\"resultdata\"><ul>";
	print "\n<li>For better visualization, residues are colored randomly in the chord diagrams. Invalid residue names like gaps (\"-\") are colored black. Long continuous black segments depict gaps in the sequence.</li>";
	print "\n<li>Short lines connecting adjacent and neighbor residues represent short-range contacts. Similarly, lines going across the plots (connecting residues with high sequence separation) represent long-range contacts.</li>";
	print "\n</ul></td></tr>";
	print "\n<tr>";
	# Results
	my $i = 1;
	foreach my $top (sort {$top_count{$a} <=> $top_count{$b}} keys %top_count){
		my $prefix = "";
		$prefix = "Random-" if $param_na_rr_fl eq "on";
		my $id = "${prefix}Top-$top";
		$id =~ s/L/L\// if $id =~ /L[0-9]/;
		print "\n<td align=\"center\">";
		print "\n<a href=\"$web_root/chord/top$top.png\" target=\"_blank\"><img src=\"$web_root/chord/top$top.png\" width=\"".(2 * $html_display_width / (scalar keys %top_count) - 5)."\" align=\"middle\"/></a>";
		if ($i == int((scalar keys %top_count)/2)){
			print "\n</td></tr><tr>";
		}
		else{
			print "\n</td>";
		}
		$i++;
		last if $i == $MAX_DIAGRAMS + 1;
	}
	print "\n</tr>";
	print "\n</table>";
	print "\n</br>";
	chdir $dir_job or confess $!;
}

sub print_roc{
	append2log("ROC computations..");
	append2log("The list of all possible contacts");
	my %all_pos_cont = (); 
	my $a = 0;
	my $b = 0;
	my $c = 0;
	my $L = length(seq_rr("native.rr"));
	for (my $i = 1; $i <= $L; $i++){
		for (my $j = $i; $j <= $L; $j++){
			next if abs($i - $j) < $min_seq_sep; 
			$a++;
			next if not defined $native_res_list{$i};
			next if not defined $native_res_list{$j};
			$b++;
			$all_pos_cont{$i." ".$j} = 0;
			# If native exists flag the contacts to 1
			if (defined $true_contacts{$i." ".$j}){
				$all_pos_cont{$i." ".$j} = 1;
				$c++;
			}
		}
	}
	append2log("All posssible contacts are $a, ignoring non-existent native residues the count is $b");
	append2log("Out of these number of contacts is $c");
	append2log("Prepare data for each source: food for ROC curves");
	# Make an empty prediction
	my %cont_list = %all_pos_cont;
	foreach (keys %cont_list){
		# when a pair is not in predicted list, it's confidence is assumed to be zero
		$cont_list{$_} = 0;
	}
	foreach my $src (sort keys %sources){
		my %this_cont_list = %cont_list;
		# update this list with available prediction confidences
		my %rows = rr_rows_ordered_in_hash("input.rr", 1000000, $src, $param_cont_type);
		foreach (keys %rows){
			my @C = split /\s+/, $rows{$_};
			$this_cont_list{$C[0]." ".$C[1]} = sprintf "%.5f", $C[4];
		}
		# now, record all contacts with true labels and available prediction confidences
		print2file("./roc/$src.txt", "PredictionConfidence\tTrueLabel\tContact");
		foreach (keys %this_cont_list){
			my @C = split /\s+/, $_;
			print2file("./roc/$src.txt", $this_cont_list{$_}."\t".$all_pos_cont{$_}."\t".$C[0]."-".$C[1]);
		}
	}
	append2log("Prepare R script to plot the ROC curves");
	print2file("roc.R", "if (!require(\"ROCR\")){");
	print2file("roc.R", "  install.packages(\"ROCR\", dependencies = TRUE)");
	print2file("roc.R", "  library(ROCR)");
	print2file("roc.R", "}");
	print2file("roc.R", "png(\"roc.png\", width=2000, height=2000, res=300)");
	print2file("roc.R", "par(mar=c(4.5, 4.5, 2, 0.2), xpd = T)");
	my $i = 1;
	foreach my $src (sort keys %sources){
		print2file("roc.R", "d$i = read.table(\"./roc/$src.txt\", sep = \"\\t\", header = T)");
		print2file("roc.R", "pr$i <- prediction(d$i\$PredictionConfidence, d$i\$TrueLabel)");
		print2file("roc.R", "pf$i <- performance(pr$i, \"tpr\", \"fpr\")");
		if ($i == 1){
			print2file("roc.R", "plot(pf$i, col = $i, lty = $i, main = \"ROC Curves\")");
		}
		else{
			print2file("roc.R", "plot(pf$i, col = $i, lty = $i, add = T)");
		}
		print2file("roc.R", "auc$i <- round(as.numeric(performance(pr$i, measure = \"auc\", x.measure = \"cutoff\")\@y.values), digits = 2)");
		$i++;
	}
	my $array_legend = "";
	my $array_color  = "";
	my $array_lty    = "";
	my $i = 1;
	foreach my $src (sort keys %sources){
		$array_legend .= "\ paste(\"$src (auc =\", auc$i,\")\")";
		$array_legend .= ", " if $i != (scalar keys %sources);
		$array_color  .= "".$i."";
		$array_color  .= ", " if $i != (scalar keys %sources);
		$array_lty    .= "".$i."";
		$array_lty    .= ", " if $i != (scalar keys %sources);
		$i++;
	}
	print2file("roc.R", "legend(\"bottomright\", c($array_legend), col = c($array_color), lty = c($array_lty))");
	append2log("Plot ROC curves - Run R..");
	system_cmd("Rscript  --no-save --no-restore roc.R");
	#system_cmd("Rscript  --no-save --no-restore --verbose roc.R 2>&1"); # for debugging

	print "\n<table class=\"results\">";
	print "\n<th colspan=13 class=\"resultsheadingrow\"><b>Section ".Roman($sectionID++)."</b> ROC Curves and AUC  <font class=\"showdescription\"><a onclick=\"toggle_it('rocaucdesc')\">see description</a></font></th>";
	# Hidden description
	print "\n<tr style=\"display:none\" id=\"rocaucdesc\"><td colspan=100 class=\"resultdata\"><ul>";
	print "\n<li>Receiver operating characteristic (ROC) curve and AUC (Area under the curve) or AUROC calculated using the R library 'ROCR'.</li>";
	print "\n</ul></td></tr>";
	print "\n<tr>";
	# Results
	print "\n<td align=\"left\">";
	my $size = $html_display_width / 2 - 10;
	print "\n<a href=\"$web_root/roc.png\" target=\"_blank\"><img src=\"$web_root/roc.png\" width=\"$size\" align=\"left\"/></a>";
	print "\n</td>";
	print "\n</tr>";
	print "\n</table>";
	print "\n</br>";
}

sub print_auc_pr{
	append2log("AUC_PR Calculations using http://mark.goadrich.com/programs/AUC/");
	chdir "$dir_job/auc_pr" or error_exit($!);
	system_cmd("cp $dir_job/roc/*.txt ./");
	foreach my $src (sort keys %sources){
		system_cmd("java -jar $APP_ROOT/auc.jar $src.txt list 2>/dev/null > auc_pr_$src.log");
		system_cmd("tail -2 auc_pr_$src.log > $src.txt");
	}
	foreach my $src (sort keys %sources){
		my $line = `head -1 $src.txt`;
		chomp $line;
		my @C = split /\s+/, $line;
		my $src_with_auc = $src."(auc=".(sprintf "%.2f", $C[9]).")";
		system_cmd("sed \"s\/\$\/\\t$src_with_auc\/\" $src.txt.pr >> pr.txt");
	}
	print2file("pr.R", "if (!require(\"ggplot2\")){");
	print2file("pr.R", "  install.packages(\"ggplot2\", dependencies = TRUE)");
	print2file("pr.R", "  library(ggplot2)");
	print2file("pr.R", "}");
	print2file("pr.R", "d1 = read.table(\"pr.txt\", sep = \"\t\", header = F)");
	print2file("pr.R", "png(\"auc_pr.png\", width=3000, height=1800, res=450)");
	print2file("pr.R", "p1 <- ggplot(data = d1, aes(d1\$V1, d1\$V2)) +");
	print2file("pr.R", "  geom_line(aes(colour = factor(d1\$V3)), lwd = 0.5) +");
	print2file("pr.R", "  xlab(\"Recall\") +");
	print2file("pr.R", "  ylab(\"Precision\") +");
	print2file("pr.R", "  theme_bw() +");
	print2file("pr.R", "  theme(legend.position = \"right\", legend.title=element_blank())");
	print2file("pr.R", "print(p1)");
	append2log("Plot AUC_PR curves - Run R..");
	system_cmd("Rscript --no-save --no-restore pr.R");
	#system_cmd("Rscript  --no-save --no-restore --verbose pr.R 2>&1"); # for debugging
	chdir $dir_job or error_exit($!);

	print "\n<table class=\"results\">";
	print "\n<th colspan=13 class=\"resultsheadingrow\"><b>Section ".Roman($sectionID++)."</b> Precision-Recall Curves and AUC-PR <font class=\"showdescription\"><a onclick=\"toggle_it('prdesc')\">see description</a></font></th>";
	# Hidden description
	print "\n<tr style=\"display:none\" id=\"prdesc\"><td colspan=100 class=\"resultdata\"><ul>";
	print "\n<li>Receiver operating characteristic (ROC) curve and AUC_PR (Area under the curve) calculated using <a href=\"http://mark.goadrich.com/programs/AUC/\" target=\"_blank\">AUCCalculator 0.2</a></li>";
	print "\n</ul></td></tr>";
	print "\n<tr>";
	# Results
	print "\n<td align=\"left\">";
	my $size = $html_display_width / 2 - 10;
	print "\n<a href=\"$web_root/auc_pr/auc_pr.png\" target=\"_blank\"><img src=\"$web_root/auc_pr/auc_pr.png\" width=\"$size\" align=\"left\"/></a>";
	print "\n</td>";
	print "\n</tr>";
	print "\n</table>";
	print "\n</br>";

}

sub print_MCC{
	append2log("MCC Calculations");
	foreach my $src (sort keys %sources){
		my %rr_selected = %{$selected_contacts{$src." "."above0.5confidence"}};
		append2log("\nWarning! There are no $param_atom_type $param_cont_type contacts with prediction confidence >0.5 for the source $src of ! Cannot proceed!") if not scalar keys %rr_selected;
		$mcc{$src} = calc_mcc(\%rr_selected);
	}
	print "\n<table class=\"results\">";
	print "\n<th colspan=13 class=\"resultsheadingrow\"><b>Section ".Roman($sectionID++)."</b> Matthew's Correlation Coefficient (MCC)  <font class=\"showdescription\"><a onclick=\"toggle_it('mccdesc')\">see description</a></font></th>";
	# Hidden description
	print "\n<tr style=\"display:none\" id=\"mccdesc\"><td colspan=100 class=\"resultdata\"><ul>";
	print "\n<li><img src=\"http://iris.rnet.missouri.edu/coneva/mcc-eqn.png\" width=\"325\"/></li>";
	print "\n<li><b>TP</b> is the count of true positives, <b>FP</b> is the count of false positives, <b>TN</b> is the count of true negatives, and <b>FN</b> is the count of false negatives.</li>";
	print "\n</ul></td></tr>";
	print "<tr><td>";
	print "\n<table>";
	print "\n<tr><td align=\"right\">RR-File</td>";
	print "\n<td></td><td align=\"right\">TP</td><td></td><td align=\"right\">FP</td>";
	print "\n<td></td><td align=\"right\">TN</td><td></td><td align=\"right\">FN</td>";
	print "\n<td></td><td align=\"right\">MCC</td></tr>";
	foreach my $src (sort keys %sources){
		print "\n<tr>";
		print "\n<td class=\"smalldataboxvalues\">$src</td>";
		my @C = split /\s+/, $mcc{$src};
		print "\n<td></td><td class=\"smalldataboxvalues\">".$C[0]."</td>";
		print "\n<td></td><td class=\"smalldataboxvalues\">".$C[1]."</td>";
		print "\n<td></td><td class=\"smalldataboxvalues\">".$C[2]."</td>";
		print "\n<td></td><td class=\"smalldataboxvalues\">".$C[3]."</td>";
		print "\n<td></td><td class=\"smalldataboxvalues\">".$C[4]."</td>";
		print "\n</tr>";
	}
	print "</table>";
	print "</td></tr>";
	print "\n</table>";
	print "\n</br>";
}

sub chk_errors_and_log_params{
	$param_rr_file  = undef if (defined $param_rr_file and length($param_rr_file) < 2);
	$param_rr_raw   = undef if (defined $param_rr_raw and length($param_rr_raw) < 2);
	$param_sec      = undef if (defined $param_sec and length($param_sec) < 2);
	$param_pdb_file = undef if (defined $param_pdb_file and length($param_pdb_file) < 2);
	$param_pdb_raw  = undef if (defined $param_pdb_raw and length($param_pdb_raw) < 2);
	$param_rr_raw   =~ s/\r//g if defined $param_rr_raw;
	# trim space from both ends
	$param_rr_raw =~ s/^\s+|\s+$//g if defined $param_rr_raw;
	$param_pdb_raw =~ s/^\s+|\s+$//g if defined $param_pdb_raw;
	# If upload file is present raw text fields should be blank
	error_exit("When uploading a RR file, the RR text field should be blank!") if defined $param_rr_raw and defined $param_rr_file;
	error_exit("When uploading a PDB file, the RR text field should be blank!") if defined $param_pdb_raw and defined $param_pdb_file;
	# Secondary structure check
	if (defined $param_sec){
		chomp $param_sec;
		for(my $i = 1; $i <= length($param_sec); $i++){
			my $char = substr $param_sec, $i-1, 1;
			if (not ($char eq "H" or $char eq "C" or $char eq "E")){
				error_exit("undefined secondary structure unit $char in $param_sec");
			}
		}
	}
	# Log all the parameters
	if (defined $param_rr_raw){
		my @L = split /\n/, $param_rr_raw;
		print2file($file_history, "rr(row1) : ".$L[0]);
	}
	print2file($file_history, "time     : ".$now);
	print2file($file_history, "rr_file  : ".$param_rr_file) if defined $param_rr_file;
	print2file($file_history, "native fl: ".$param_na_rr_fl) if defined $param_na_rr_fl;
	print2file($file_history, "pdb(row1): ".substr($param_pdb_raw, 0, 100)) if defined $param_pdb_raw;
	print2file($file_history, "pdb_file : ".$param_pdb_file) if defined $param_pdb_file;
	print2file($file_history, "sec      : ".$param_sec) if defined $param_sec;
	# Set distance threshold to undefined because it must be initialized from input RR files
	if($param_na_rr_fl ne "on"){
		$d_threshold = undef;
	}
}

sub create_files{
	append2log("Creating files");
	# Param 2 : Raw RR file
	if (defined $param_rr_raw){
		if($param_rr_raw !~ /^http/){
			append2log("Content in Text area");
			print2file("input_original.rr", $param_rr_raw);
			$flag_1_src = 1;
			# Select only a certain subset of contacts, because usually too many contacts are predicted!
			$rr_sequence = seq_rr("input_original.rr") if not defined $rr_sequence;
			my %rr_rows = rr_rows_ordered_in_hash("input_original.rr", ($MAX_CONTACTS * length($rr_sequence)), "all", $param_cont_type);
			print2file("input.rr", $rr_sequence);
			foreach my $row (sort {$a <=> $b} keys %rr_rows){
				print2file("input.rr", $rr_rows{$row});
			}
			sort_rr_file_by_confidence("$dir_job/input.rr");
			$input_rr = "$dir_job/input.rr";
		}
		else{
			# URL of CASP10 RR location can be accepted
			system_cmd("wget $param_rr_raw");
			my @suffixlist = qw/.tgz .tar.gz .zip/;
			my($filename, $dirs, $suffix) = fileparse($param_rr_raw);
			my $basename_zip = basename($filename, @suffixlist);
			if($filename =~ /zip/){
				system_cmd("unzip $filename -d ./ >> $file_log");
				make_rr_from_input_folder($basename_zip);
			}
			elsif($filename =~ /tar/ or $filename =~ /tgz/){
				system_cmd("tar zxvf $filename >> $file_log");
				make_rr_from_input_folder($basename_zip);
			}
			elsif($filename =~ /RR/ and $filename =~ /_1/){
				sort_rr_file_by_confidence($filename);
				$flag_1_src = 1;
				$rr_sequence = seq_rr($filename);
				# Check if the files contain contact types of interest
				# This also fixes the error of R1 > R2;
				my %rr_rows = rr_rows_ordered_in_hash($filename, ($MAX_CONTACTS * length($rr_sequence)), "all", $param_cont_type);
				error_exit("Input RR file does not contain any $cont_type_desc contacts") if (scalar keys %rr_rows) < 1;
				print2file("input.rr", $rr_sequence);
				foreach my $row (sort {$a <=> $b} keys %rr_rows){
					print2file("input.rr", $rr_rows{$row});
				}
				$input_rr = "$dir_job/input.rr";
			}
			else{
				system_cmd("cp $filename input.rr");
				$input_rr = "$dir_job/input.rr";
				$rr_sequence = seq_rr($input_rr);
				sort_rr_file_by_confidence($input_rr);
				$flag_1_src = 1;
				# Check if the files contain contact types of interest
				# This also fixes the error of R1 > R2;
				my %rr = rr_rows_ordered_in_hash($input_rr, 100000, "all", $param_cont_type);
				error_exit("Input RR file does not contain any $cont_type_desc contacts") if (scalar keys %rr) < 1;
				system_cmd("rm -f temp.rr");
				print2file("temp.rr", seq_rr($input_rr));
				foreach (sort {$a <=> $b} keys %rr){
					print2file("temp.rr", $rr{$_});
				}
				system_cmd("mv temp.rr $input_rr");
			}
		}
		$input_rr = "$dir_job/input.rr";
		# Some files have leading white spaces
		system_cmd("sed -i 's/^ *//' $input_rr");
	}
	# Param 3: Upload an RR file
	if (defined $param_rr_file){
		error_exit("Filename contains invalid characters") if $param_rr_file !~ m/[$safe_filename_chars]+/;
		# copy the file as it is
		append2log("Received input uploaded file $param_rr_file");
		open INPUT, ">$param_rr_file" or confess $!;
		while(<$param_rr_file>) {
			print INPUT $_;
		}
		close INPUT;
		# check if it is a text or a zipped file
		my @suffixlist = qw/.tgz .tar.gz .zip/;
		my($filename, $dirs, $suffix) = fileparse($param_rr_file);
		my $basename_zip = basename($filename, @suffixlist);
		if($filename =~ /zip/){
			system_cmd("unzip $filename -d ./ >> $file_log");
			make_rr_from_input_folder($basename_zip);
		}
		elsif($filename =~ /tar/ or $filename =~ /tgz/){
			system_cmd("tar zxvf $filename >> $file_log");
			make_rr_from_input_folder($basename_zip);
		}
		else{
			system_cmd("cp $param_rr_file input.rr");
			$input_rr = "$dir_job/input.rr";
			$rr_sequence = seq_rr($input_rr);
			$flag_1_src = 1;
			sort_rr_file_by_confidence($input_rr);
		}
	}
	# Param 4: Raw PDB file
	if (defined $param_pdb_raw){
		# URL of pdb file 
		if($param_pdb_raw =~ /^http/){
			system_cmd("wget $param_pdb_raw");
			my($filename, $dirs, $suffix) = fileparse($param_pdb_raw);
			my $file = basename($filename); 
			error_exit("ERROR! Downloaded pdb file $file does not exist!") if not -f $file;
			error_exit("ERROR! Downloaded pdb file $file is empty!") if not -s $file;
			system_cmd("mv $file native.pdb");
		}
		# PDB id supplied
		elsif($param_pdb_raw =~ /^id/){
			my @C = split /:/, lc($param_pdb_raw);
			my $chain_id;
			if (length($C[1]) > 4){
				$chain_id = uc(substr $C[1], 4, 1);
				$C[1] = substr $C[1], 0, 4; 
			}
			$param_pdb_raw = "http://www.rcsb.org/pdb/files/".$C[1].".pdb";
			system_cmd("wget $param_pdb_raw");
			my($filename, $dirs, $suffix) = fileparse($param_pdb_raw);
			my $file = basename($filename); 
			if ((not -f $file) or  (not -f $file)){
				error_exit("ERROR! Attempted downloading $param_pdb_raw</br> Looks like it does not exist!</br>Could not find $file in downloads location!");
			}
			if($chain_id){
				system("grep '^ATOM.\\{17\\}$chain_id' $file > native.pdb");
				if(not -s "native.pdb"){
					error_exit("ERROR! Chain $chain_id does not exist in $param_pdb_raw!");
				}
				# for summary display it is better to add chain information
				$param_pdb_raw = $param_pdb_raw." (chain $chain_id)";
			}
			else{
				system_cmd("mv $file native.pdb");
			}
		}
		else{
			# Content in textarea
			print2file("native.pdb", $param_pdb_raw);
		}
		$native_pdb = "native.pdb";
	}
	# Param 5: Upload a PDB file
	if (defined $param_pdb_file){
		append2log("Creating native PDB from uploaded link contents");
		error_exit("Filename contains invalid characters") if $param_pdb_file !~ m/[$safe_filename_chars]+/;
		open PDB, ">native.pdb" or error_exit ($!);
		while(<$param_pdb_file>) {
			print PDB $_;
		}
		close PDB;
		$native_pdb = "native.pdb";
	}
	# Check the input native PDB
	if (-f $native_pdb){
		my $atom_row_count = 0;
		open CHAIN, $native_pdb or confess $!;
		while(<CHAIN>){
			next if $_ !~ m/^ATOM/;
			$atom_row_count++;
		}
		close CHAIN;
		error_exit("It looks like the native PDB file has too few atom rows! Please check your input PDB file!") if $atom_row_count < 10;
	}
	# Param 6: Secondary structure
	if (defined $param_sec){
		print2file("input.ss", ">$remote_addr");
		print2file("input.ss", $param_sec);
		$input_sec = "input.ss";
	}
	# Param 7: Native RR analysis
	if($param_na_rr_fl eq "on"){
		#reindex_chain($native_pdb, 1, "reindexed.pdb");
		#$native_pdb = "reindexed.pdb";
		pdb2rr($native_pdb, "input.rr", $min_seq_sep, $d_threshold, $param_atom_type);
		$input_rr = "$dir_job/input.rr";
		$rr_sequence = seq_rr($input_rr);
		append2log("Native RR flag is on, using sequence $rr_sequence");
	}
}

sub make_rr_from_input_folder{
	my $basename_zip = shift;
	append2log("Making input.rr from foler $basename_zip");
	system_cmd("rename ' ' '_' ./$basename_zip/*");
	system_cmd("rename ' ' '_' ./$basename_zip/*");
	system_cmd("rename ' ' '_' ./$basename_zip/*");
	system_cmd("rename ' ' '_' ./$basename_zip/*");
	system_cmd("rename ' ' '_' ./$basename_zip/*");
	system_cmd("rename ' ' '_' ./$basename_zip/*");
	my @rr_files = <./$basename_zip/*>;
	error_exit("Could not find any RR files in the <b>$basename_zip</b> folder after unzipping!</br>Please put the files into a folder and zip the folder!</br>Also please make sure that the zipped file name is same as the name of the folder that is being zipped!") if not -f $rr_files[0];
	if (rr_has_seq($rr_files[0]) and defined $rr_sequence){
		error_exit("RR file must have some sequence input !") if not defined seq_rr($rr_files[0]);
		error_exit("RR input sequence and file sequence do not match!") if seq_rr($rr_files[0]) ne $rr_sequence;
	}
	$rr_sequence = seq_rr($rr_files[0]) if not defined $rr_sequence; 
	print2file("input.rr", $rr_sequence);
	$flag_1_src = 0 if defined $rr_files[1];
	foreach my $file(@rr_files) {
		# Check that all input sequences are same, except for the case 3 when native sequence is combined with input sequence
		if (!($file =~ /input\/native/ or $rr_files[0] =~ /input\/native/)){
			if ($rr_sequence ne seq_rr($file)){
				error_exit("Sequence mismatch error!</br>Sequence of ".$rr_files[0]." and $file do not match!<br>$rr_sequence [".$rr_files[0]."]</br>".seq_rr($file)."[$file]");
			}
		}
		my $base_rr = basename($file);
		sort_rr_file_by_confidence($file);
		append2log("Reading $param_cont_type input from $base_rr");
		my %rr_rows = rr_rows_ordered_in_hash($file, ($MAX_CONTACTS * length($rr_sequence)), "all", $param_cont_type);
		append2log("Appending contacts from $base_rr");
		foreach my $row (sort {$a <=> $b} keys %rr_rows){
			my @C = split /\s+/, $rr_rows{$row};
			error_exit("Contact row ".$rr_rows{$row}." has errors!") if not defined $C[4];
			print2file("input.rr", $rr_rows{$row}) if $flag_1_src == 1;
			# Add a sixth column in case of multiple input files
			print2file("input.rr", $rr_rows{$row}." $base_rr") if $flag_1_src == 0;
		}
	}
	$input_rr = "$dir_job/input.rr";
}

sub chk_errors_input_files{
	# Sequence related validations
	if (defined $rr_sequence){
		$rr_sequence = seq_fasta($fasta) if (-f $fasta and not defined $rr_sequence);
		error_exit("Sequence mismatch between Input Sequence and RR!!") if $rr_sequence ne seq_rr($input_rr);
		if (-f $native_pdb){
			system_cmd("cp /var/www/html/coneva/header.php sequence_compare.html");
			print2file("sequence_compare.html", "<table align=\"center\">");
			print2file("sequence_compare.html", "<tr><td class=\"coordinationnumber\">");
			if ($rr_sequence ne seq_chain_with_gaps($native_pdb)){
				print2file("sequence_compare.html", "<b><font color=\"red\">WARNING!</font> Sequence mismatch between Input Sequence and PDB! Some predicted contacts not found in PDB may be igonored!</b></br>");
			}
			else{
				print2file("sequence_compare.html", "<b>Sequence comparison of Input Sequence and PDB!</b></br>");
			}
			print2file("sequence_compare.html", "$rr_sequence [Input RR Sequence]</br>");
			print2file("sequence_compare.html", "".seq_chain_with_gaps($native_pdb)." [PDB Sequence]");
			print2file("sequence_compare.html", "</td></tr></table>");
			system_cmd("cat /var/www/html/coneva/footer.php >> sequence_compare.html");
		}
		if (-f $input_sec){
			error_exit("Sequence Length mismatch between Input Sequence and SS!") if length($rr_sequence) != length(seq_fasta($input_sec));
		}
	}
=notneeded
	if (not defined $rr_sequence){
		$rr_sequence = seq_rr($input_rr);
		if (-f $native_pdb){
			if ($rr_sequence ne seq_chain_with_gaps($native_pdb)){
				system_cmd("cp /var/www/html/coneva/header.php sequence_compare.html");
				print2file("sequence_compare.html", "<table>");
				print2file("sequence_compare.html", "<tr><td class=\"coordinationnumber\">");
				print2file("sequence_compare.html", "<b><font color=\"red\">WARNING!</font> Sequence mismatch between Input Sequence and PDB!</b></br>");
				print2file("sequence_compare.html", "$rr_sequence [Input RR Sequence]</br>");
				print2file("sequence_compare.html", "".seq_chain_with_gaps($native_pdb)." [PDB Sequence]");
				print2file("sequence_compare.html", "</td></tr></table>");
				system_cmd("cat /var/www/html/coneva/footer.php >> sequence_compare.html");
			}
		}
		if (-f $input_sec){
			error_exit("Sequence Length mismatch between Input Sequence and SS!") if length($rr_sequence) != length(seq_fasta($input_sec));
		}
	}
=cut
	error_exit("Cannot find sequence anywhere!!") if not defined $rr_sequence;
	chk_errors_seq($rr_sequence);
	# Add sequence to the top of RR file, if one does not exist
	if (not defined seq_rr($input_rr)){
		system_cmd("rm -f temp.rr");
		print2file("temp.rr", $rr_sequence);
		system_cmd("cat $input_rr >> temp.rr");
		system_cmd("mv temp.rr $input_rr");
	}
	error_exit("No input RR file!") if not -f $input_rr;
	
	#checking rr errors
	open RR, $input_rr or error_exit($!);
	my $no_rows = 1;
	while(<RR>){
		my $line = $_;
		next if length($line) < 2;
		next if $line =~ /^[A-Z]/;
		next if $line =~ /^-/;
		next if $line =~ /^>/;
		my @C = split /\s+/, $line;
		error_exit("Column 1 not defined in [row: $line] in rr input") if not defined $C[0];
		error_exit("Column 2 not defined in [row: $line] in rr input") if not defined $C[1];
		error_exit("Column 3 not defined in [row: $line] in rr input") if not defined $C[2];
		error_exit("Column 4 not defined in [row: $line] in rr input") if not defined $C[3];
		error_exit("Column 5 not defined in [row: $line] in rr input") if not defined $C[4];
		error_exit("First residue in [row: $line] is smaller than the sequence length! Something is wrong! Please check the input RR file!") if $C[0] > length($rr_sequence);
		error_exit("Second residue in [row: $line] is smaller than the sequence length! Something is wrong! Please check the input RR file!") if $C[1] > length($rr_sequence);
		$no_rows = 0;
	}
	close RR;
	error_exit("Input RR file has 0 rows that satisfy the input contact definition! Please (i) check the input RR file (ii) try loosening the contact definition like decreasing the minimum sequence separation or increase the distance threshold!") if $no_rows;
	# Find the distance threshold
	open RR, $input_rr or error_exit($!);
	while(<RR>){
		my $line = $_;
		next if length($line) < 2;
		next if $line =~ /^[A-Z]/;
		next if $line =~ /^-/;
		next if $line =~ /^>/;
		my @C = split /\s+/, $line;
		$d_threshold = $C[3];
		last;
	}
	close RR;
}

sub print_footer{
	open FOOTER, "/var/www/html/coneva/footer.php"; 
	while(<FOOTER>){
		print $_;
	}
	close FOOTER;
	exit;
}

sub error_exit{
	my $message = shift;
	print2file($file_history, $message);
	print "\n<table width=\"1000\" align=\"center\" cellspacing=\"0\" cellpadding=\"5\">";
	print "\n<tr bgcolor=\"#CCCCDD\"><td><b><font color=\"red\">ERROR!!</font></b></td><td></td></tr>";
	print "\n<tr><td>";
	print "\n$message";
	print "\n</td></tr></table>";
	print_footer();
	exit;
}

sub aa_color{
	my $c = shift;
	my %aacol = ();
	$aacol{"A"} = "#e00000";
	$aacol{"N"} = "#c00000";
	$aacol{"C"} = "#f00030";
	$aacol{"Q"} = "#f00000";
	$aacol{"H"} = "#c00000";
	$aacol{"L"} = "#e00060";
	$aacol{"M"} = "#c00010";
	$aacol{"P"} = "#800000";
	$aacol{"T"} = "#800030";
	$aacol{"Y"} = "#700010";
	$aacol{"R"} = "#400000";
	$aacol{"D"} = "#7000c0";
	$aacol{"E"} = "#3000a0";
	$aacol{"G"} = "#200080";
	$aacol{"I"} = "#a000a0";
	$aacol{"K"} = "#700070";
	$aacol{"F"} = "#500060";
	$aacol{"S"} = "#e00020";
	$aacol{"W"} = "#c00000";
	$aacol{"V"} = "#a00000";
	$aacol{"X"} = "#500050";
	return $aacol{$c};
}

sub rr_rows_ordered_in_hash{
	my $file_rr = shift;
	my $count = shift;
	my $source = shift;
	my $type = shift;
	my $confidence_thres = shift;
	confess "Input file not defined" if not defined $file_rr;
	confess "Input file $file_rr does not exist!" if not -f $file_rr;
	confess "No File RR!" if not -f $file_rr;
	confess "No contact count!" if not defined $count;
	confess "No contact source!" if not defined $source;
	confess "No contact type" if not defined $type;
	confess "Invalid type" if not ($type eq "everything" or $type eq "all" or $type eq "short" or $type eq "medium" or $type eq "long");
	my %rr = ();
	my $i = 1;
	my %i_for_each_source = ();
	# Find all Contact sources
	if ($source eq "all"){
		open RR, $file_rr or confess $!;
		while (<RR>){
			next unless $_ =~ /[0-9]/;
			chomp $_;
			$_ =~ s/\r//g;
			$_ =~ s/^\s+//;
			next unless $_ =~ /^[0-9]/;
			my @C = split /\s+/, $_ ;
			last if not defined $C[5];
			last if (scalar keys %i_for_each_source) >= $MAX_SOURCES;
			$i_for_each_source{$C[5]} = 1;
		}
		close RR;
	}
	open RR, $file_rr or error_exit($!);
	while(<RR>){
		my $row = $_;
		next unless $row =~ /[0-9]/;
		chomp $row;
		$row =~ s/\r//g;
		$row =~ s/^\s+//;
		next unless $row =~ /^[0-9]/;
		my @C = split /\s+/, $row ;
		error_exit("Expecting a pair in row [".$row."]!\n") if (not defined $C[0] || not defined $C[1]);
		error_exit("Confidence column not defined in row [".$row."] in file <b>$file_rr</b>! </br>Please make sure that the input RR file is in 5-column format!") if not defined $C[4];
		# Fix order 
		if ($C[0] > $C[1]){
			$row = $C[1]." ".$C[0]." ".$C[2]." ".$C[3]." ".$C[4];
			$row = $C[1]." ".$C[0]." ".$C[2]." ".$C[3]." ".$C[4]." ".$C[5] if defined $C[5];
		}
		if (defined $confidence_thres){
			next if $C[4] < $confidence_thres;
		}
		# Select only LR, MR, SR or all
		my $d = abs($C[0]-$C[1]);
		if ($type eq "long"){
			next if $d < $long_range_min;
		}
		if ($type eq "medium"){
			next if $d < $medum_range_min;
			next if $d > $medium_range_max;
		}
		if ($type eq "short"){
			next if $d < $short_range_min;
			next if $d > $short_range_max;
		}
		next if ($d < $short_range_min and $type eq "all");
		# Select only those matching the source
		$C[5] = "all" if not defined $C[5];
		if ($source eq "all" and $C[5] eq "all"){
			$rr{$i} = $row;
			$i++;
			last if $i > $count;
		}
		elsif($source ne "all" and $C[5] eq "all"){
			confess "Specified some specific source but there is a whose source is not specified";
		}
		elsif($source eq "all" and $C[5] ne "all"){
			next if not defined $i_for_each_source{$C[5]};
			next if $i_for_each_source{$C[5]} > $count;
			$i++;
			$rr{$i} = $row;
			$i_for_each_source{$C[5]}++;
		}
		else{
			next if $C[5] ne $source;
			$rr{$i} = $row;
			$i++;
			last if $i > $count;
		}
	}
	close RR;
	append2log("Subroutine: rr_rows_ordered_in_hash :: File: $file_rr Type: $type Count: $count Source: $source -> Returning ".(scalar keys %rr)." rows");
	return %rr;
}

sub sort_rr_file_by_confidence{
	my $file_rr = shift;
	error_exit("No File RR!") if not -f $file_rr;
	my $seq = undef;
	if(rr_has_seq($file_rr)){
		$seq = seq_rr($file_rr);
	}
	system_cmd("rm -f sorted.rr");
	# Keep contact rows only
	system_cmd("sed -i '/^[A-Z]/d' $file_rr");
	system_cmd("sed -i '/^-]/d' $file_rr");
	# Some files have leading white spaces
	system_cmd("sed -i 's/^ *//' $file_rr");
	# Stable sort with -s option, i.e. maintain order in case confidence are equal
	# Also using -g option instead of -n because some predictors have exponential values in confidence
	system_cmd("sort -gr -s -k5 $file_rr > sorted.rr");
	system_cmd("rm -f $file_rr");
	print2file($file_rr, $seq) if defined $seq;
	system_cmd("cat sorted.rr >> $file_rr");
	system_cmd("rm sorted.rr");
}

sub chk_rr_not_in_native{
	if(-f $native_pdb){
		my $ignore_count = 0;
		my $total = 0;
		system_cmd("rm -f temp.rr");
		print2file("temp.rr", $rr_sequence);
		my $flag_has_at_least1 = 0;
		open RR, $input_rr or error_exit($!);
		while(<RR>){
			my $line = $_;
			chomp $line;
			next if $line =~ /^[A-Z]/;
			next if $line =~ /^-/;
			$total++;
			my @R = split /\s+/, $line;
			if (not defined $native_res_list{$R[0]} or not defined $native_res_list{$R[1]}){
				$ignore_count++;
				print2file("./not-in-native/not_in_native_".$R[5].".rr", $line) if defined $R[5];
				print2file("./not-in-native/not_in_native_all.rr", $line) if not defined $R[5];
				next;
			}
			print2file("temp.rr", $line);
			$flag_has_at_least1 = 1;
		}
		close RR;
		system_cmd("mv $input_rr input_unfiltered.rr");
		if (!$flag_has_at_least1){
			error_exit("None of the input predicted contacts were found in the native pdb! Please check the input native pdb file and the RR file. Some reindexing is probably needed!");
		}
		system_cmd("mv temp.rr $input_rr");
		if ($ignore_count){
			$total_input_rr = $total;
			$ignored_rr_cnt = $ignore_count;
		}
	}
}

sub chk_errors_seq{
	my $rr_sequence = shift;
	if (length($rr_sequence) < $MIN_SEQ_LIMIT){
		error_exit("Sequence $rr_sequence is undefined or too short!");
	}
	if (length($rr_sequence) > $MAX_SEQ_LIMIT){
		error_exit("Sequence $rr_sequence is too long.</br>$MAX_SEQ_LIMIT is the limit!");
	}
	for(my $i = 1; $i <= length($rr_sequence); $i++){
		my $char = substr $rr_sequence, $i-1, 1;
		if (not defined $AA1TO3{$char}){
			error_exit("undefined amino acid $char in $rr_sequence");
		}
	}
}

sub calc_jaccard_similarity{
	my $rr1name = shift;
	my $rr2name = shift;
	my $rr_set1 = shift;
	my $rr_set2 = shift;
	my %rr1input = %{$rr_set1};
	my %rr2input = %{$rr_set2};
	my %rr1 = ();
	my %rr2 = ();
	my %union = ();
	my %intersect = ();
	foreach (keys %rr1input){
		my @C = split /\s+/, $rr1input{$_};
		$rr1{$C[0]." ".$C[1]} = 1;
	}
	foreach (keys %rr2input){
		my @C = split /\s+/, $rr2input{$_};
		$rr2{$C[0]." ".$C[1]} = 1;
	}
	# add to union set
	foreach (keys %rr1){
		my @C = split /\s+/, $_;
		my $flag_exists = 0;
		for (my $i = -($JN); $i <= $JN; $i++){
			for (my $j = -($JN); $j <= $JN; $j++){
				if (defined $union{($C[0] + $i)." ".($C[1] + $j)}){
					$flag_exists = 1;
					last;
				}
			}
		}
		$union{$_} = 1 if ! $flag_exists;
	}
	foreach (keys %rr2){
		my @C = split /\s+/, $_;
		my $flag_exists = 0;
		for (my $i = -($JN); $i <= $JN; $i++){
			for (my $j = -($JN); $j <= $JN; $j++){
				if (defined $union{($C[0] + $i)." ".($C[1] + $j)}){
					$flag_exists = 1;
					last;
				}
			}
		}
		$union{$_} = 1 if ! $flag_exists;
	}
	# add to intersection set
	foreach (keys %union){
		my @C = split /\s+/, $_;
		my $flag_exists_rr1 = 0;
		my $flag_exists_rr2 = 0;
		for (my $i = -($JN); $i <= $JN; $i++){
			for (my $j = -($JN); $j <= $JN; $j++){
				$flag_exists_rr1 = 1 if defined $rr1{($C[0] + $i)." ".($C[1] + $j)};
				$flag_exists_rr2 = 1 if defined $rr2{($C[0] + $i)." ".($C[1] + $j)};
				last if ($flag_exists_rr1 and $flag_exists_rr2);
			}
		}
		$intersect{$_} = 1 if ($flag_exists_rr1 and $flag_exists_rr2);
	}
	my $js = sprintf "%.2f", (scalar keys %intersect)/((scalar keys %union));
	append2log("jaccard similarity between $rr1name and $rr2name: #set1 = ".(scalar keys %rr1)." #set2 = ".(scalar keys %rr2)." #Union = ".(scalar keys %union)." #Intersect = ".(scalar keys %intersect)." #JaccardSimilariy = $js");
	return $js;
}

sub calc_mcc{
	my $rrhash = shift;
	my %rr = %{$rrhash};
	if (not scalar keys %rr){
		return "- - - - -";
	}
	my %distances = rrhash2dist(\%rr, \%min_true_dist);
	error_exit("Distances could not be calculated for selected contacts! Something went wrong!") if not scalar keys %distances;
	my $satisfied = 0;
	foreach (sort {$distances{$a} <=> $distances{$b}}keys %distances){
		my @R = split /\s+/, $_;
		$satisfied++ if $distances{$_} <= $R[3];
	}
	my $tp = $satisfied;
	my $fp = (scalar keys %rr) - $satisfied;
	my $tn = (scalar keys %min_true_dist) - (scalar keys %true_contacts) - $fp;
	my $fn = (scalar keys %true_contacts) - $satisfied;
	my $mcc = ($tp*$tn - $fp*$fn)/sqrt(($tp+$fp)*($tp+$fn)*($tn+$fp)*($tn+$fn)); 
	# For debugging purposes
	#return sprintf "%.2f %.2f(".(scalar keys %rr)."-$satisfied) \
	#%.2f(".(scalar keys %min_true_dist)."-".(scalar keys %true_contacts)."-".$fp.") %.2f(".(scalar keys %true_contacts)."-$satisfied) %.2f", $tp, $fp, $tn, $fn, $mcc;
	return sprintf "%d %d %d %d %.2f", $tp, $fp, $tn, $fn, $mcc;
}

sub calc_precision{
	my $rrhash = shift;
	my %rr = %{$rrhash};
	error_exit("Cannot calculate precision! Empty selected contacts") if not scalar keys %rr;
	my %distances = rrhash2dist(\%rr, \%min_true_dist);
	error_exit("Distances could not be calculated for selected contacts! Something went wrong!") if not scalar keys %distances;
	my $satisfied = 0;
	foreach (sort {$distances{$a} <=> $distances{$b}}keys %distances){
		my @R = split /\s+/, $_;
		$satisfied++ if $distances{$_} <= $R[3];
	}
	return sprintf "%.2f", 100 * ($satisfied/(scalar keys %distances));
}

sub calc_coverage{
	my $rrhash = shift;
	my %rr = %{$rrhash};
	error_exit("Cannot calculate precision! Empty selected contacts") if not scalar keys %rr;
	my %rrpdb = rr_rows_ordered_in_hash($native_rr, 100000, "all", $param_cont_type);
	my $total = scalar keys %rrpdb;
	my $covered = 0;
	my %two_col_rr = ();
	foreach (keys %rr){
		my @R = split /\s+/, $rr{$_};
		$two_col_rr{$R[0]." ".$R[1]} = 1;
	}
	foreach (sort keys %rrpdb){
		my @C = split /\s+/, $rrpdb{$_};
		my $pair = $C[0]." ".$C[1]; 
		$covered++ if defined $two_col_rr{$pair};
	}
	append2log("Covered: $covered of Total: $total");
	return sprintf "%.2f", 100 * ($covered/$total);
}

sub calc_xd{
	my $rrhash = shift;
	my $xd_calc = 0;
	my %rr = %{$rrhash};
	my %pred_cont_dist = rrhash2dist(\%rr, \%min_true_dist);
	confess "ERROR! Empty input rr!" if not scalar keys %pred_cont_dist;
	for(my $i = 1; $i <= 15; $i++){
		my $dmin = 4 * ($i - 1);
		my $dmax = 4 * $i;
		my $ppi = count_in_range_xd($dmin, $dmax, \%pred_cont_dist)/(scalar keys %pred_cont_dist);
		my $pai = fraction_of_min_true_dists($dmin, $dmax);
#		my $sum = 100 * ($ppi - $pai)/(0.25 * $dmax);
		my $sum = 100 * ($ppi - $pai)/($i);
		$xd_calc += $sum;
	#	append2log(sprintf "$i $dmin $dmax %d %.3f %.3f %.3f %.3f", count_in_range_xd($dmin, $dmax, \%pred_cont_dist), $ppi, $pai, $sum, $xd_calc);
	}
	return sprintf "%.2f", $xd_calc;
}

sub rrhash2dist{
	# If input is empty output will be empty
	my $rrhash = shift;
	my %rr_hash = %{$rrhash};
	my %output = ();
	foreach (keys %rr_hash){
		my @C = split /\s+/, $rr_hash{$_};
		next if not defined $min_true_dist{$C[0]." ".$C[1]};
		$output{$rr_hash{$_}} = $min_true_dist{$C[0]." ".$C[1]};
	}
	return %output;
}

sub calc_spread{
	# Spread: distribution (spread) of the contacts along the chain and 
	# over the structure of the protein, by measuring the mean of the distance
	# from every experimental (crystal structure) contact to the nearest predicted contact 
	## We compute the Euclidean 2D distance between the contact map of the corresponding 
	## crystal structure, and the contact map consisting of the top N-scoring residue pairs. 
	## For each residue pair, separated by more than five residues in sequence, we compute the 
	## distance to the nearest high-scoring residue pair (for instance, the nearest ‘red star’ 
	## in the contact map, in the case of EIC pairs). For each set of N-scoring residue pairs we 
	## calculate the mean of the distances for all contacts to the nearest contact in the crystal 
	## structure. Plotted is the mean spread for each Nc for 4 methods, across all 15 proteins.
	my $file_rr = shift;
	my %predicted = rrfile_to_r1r2hash($file_rr, 1, 100000);
	confess "ERROR! Empty input rr!" if not scalar keys %predicted;
	my $spread = 0;
	confess "There are no contacts from the native!" if not scalar keys %true_contacts;
	foreach my $nat(keys %true_contacts){
		my $min = 100000;
		foreach my $pred(keys %predicted){
			my @N = split /\s+/, $nat;
			my @P = split /\s+/, $pred;
			my $d = sqrt(($N[0]-$P[0])*($N[0]-$P[0])+($N[1]-$P[1])*($N[1]-$P[1]));
			$min = $d if $d < $min;
		}
		$spread += $min;
	}
	return sprintf "%.2f", $spread/(scalar keys %true_contacts);
}

sub calc_fp{
	# In reference to the EVFOLD paper
	my $rrhash = shift;
	my %rr = %{$rrhash};
	error_exit("Cannot calculate fp! Empty selected contacts") if not scalar keys %rr;
	my %distances = rrhash2dist(\%rr, \%min_true_dist);
	error_exit("Distances could not be calculated for selected contacts! Something went wrong!") if not scalar keys %distances;
	confess "\nERROR! Native must have been loaded already!" if not scalar keys %min_true_dist;
	my $deviation = 0;
	foreach (keys %rr){
		my @P = split /\s+/, $rr{$_};
		confess "One of the residues in ".$P[0]." ".$P[1]." is not defined in native_pdb!" if not defined $min_true_dist{$P[0]." ".$P[1]};
		my $d = $min_true_dist{$P[0]." ".$P[1]};
		if ($d < $P[3]){
			$deviation += 0.0;
		}
		else{
			$deviation += ($d - $P[3]);
		}
	}
	return sprintf "%.2f", $deviation/(scalar keys %rr);
}

sub count_in_range_xd{
	my $dmin = shift;
	my $dmax = shift;
	my $ref_hash = shift;
	my %cont_dist = %{$ref_hash};
	my $count = 0;
	foreach my $c (keys %cont_dist){
		$count++ if ($cont_dist{$c} >= $dmin and $cont_dist{$c} <= $dmax);
	}
	return $count;
}

sub fraction_of_min_true_dists{
	my $dmin = shift;
	my $dmax = shift;
	my $count = 0;
	foreach (keys %min_true_dist){
		# Keeping < $dmax to match CASP numbers
		$count++ if ($min_true_dist{$_} >= $dmin and $min_true_dist{$_} < $dmax);
	}
	return $count/(scalar keys %min_true_dist);
}

sub load_sources{
	my $rr = shift;
	# Contact sources
	my %rr = rr_rows_ordered_in_hash($rr, 1000000, "all", "all");
	%sources = ();
	foreach (sort {$a <=> $b} keys %rr){
		my @C = split /\s+/, $rr{$_};
		last if not defined $C[5];
		$sources{$C[5]} = 1;
		last if (scalar keys %sources) >= $MAX_SOURCES;
	}
	if (not scalar keys %sources){
		$sources{"all"} = 1;
	}
	# load contacts into memory
	foreach my $top (sort {$top_count{$a} <=> $top_count{$b}} keys %top_count){
		foreach my $src (sort keys %sources){
			my %rr_selected = rr_rows_ordered_in_hash($input_rr, $top_count{$top}, $src, $param_cont_type);
			$selected_contacts{$src." ".$top} = \%rr_selected;
			my %abc = %{$selected_contacts{$src." ".$top}};
			error_exit("No contacts were found for file/type [$src]! File cannot be empty!") if not scalar keys %abc;
		}
	}
	foreach my $src (sort keys %sources){
		my %rr_selected = rr_rows_ordered_in_hash($input_rr, 100000, $src, $param_cont_type, 0.5);
		$selected_contacts{$src." "."above0.5confidence"} = \%rr_selected;
	}
	# For easier plotting..
	foreach my $s1 (sort keys %sources){
		my $short1 = substr $s1, 0, 12;
		foreach my $s2 (sort keys %sources){
			next if $s1 eq $s2;
			my $short2 = substr $s2, 0, 12;
			if($short1 eq $short2){
				error_exit("ERROR! Two rr files $s1 and $s2 have same first-12 characters ($short1 and $short2)! Please rename them so that they have unique first-12 characters!");
			}
		}
	}
}

sub load_top_count{
	my $rr = shift;
	# Top Contact subsets for analysis
	my $sequence = seq_rr($rr);
	# CASP uses Native's Length as reference
	$sequence = seq_chain($native_pdb) if -f $native_pdb;
	$top_count{"5"}  = 5;
	$top_count{"L10"} = int(0.1 * length($sequence) + 0.5);
	$top_count{"L5"}  = int(0.2 * length($sequence) + 0.5);
	$top_count{"L2"}  = int(0.5 * length($sequence) + 0.5);
	$top_count{"L"}   = length($sequence);
	$top_count{"2L"}  = int(2.0 * length($sequence));
	$top_count_order{"5"}    = 1;
	$top_count_order{"L10"}  = 2;
	$top_count_order{"L5"}   = 3;
	$top_count_order{"L2"}   = 4;
	$top_count_order{"L"}    = 5;
	$top_count_order{"2L"}   = 6;
}

sub print_contact_maps{
	append2log("Plot maps");
	chdir "$dir_job/cmap" or confess $!;
	foreach my $top (sort {$top_count{$a} <=> $top_count{$b}} keys %top_count){
		my $seq = undef;
		$seq = seq_rr($input_rr);
		system_cmd("rm -f temp.rr");
		print2file("temp.rr", $seq);
		if(-f "$dir_job/$native_pdb"){
			my %rrn = rr_rows_ordered_in_hash("$dir_job/$native_rr", 100000, "all", "all");
			foreach my $row (sort {$a <=> $b} keys %rrn){
				print2file("temp.rr", $rrn{$row}." PDB");
			}
		}
		my %rr = rr_rows_ordered_in_hash($input_rr, $top_count{$top}, "all", $param_cont_type);
		foreach my $row (sort {$a <=> $b} keys %rr){
			my @C = split /\s+/, $rr{$row};
			$rr{$row} .= " input" if not defined $C[5];
			print2file("temp.rr", $rr{$row});
		}
		my $prefix = "";
		$prefix = "Random-" if $param_na_rr_fl eq "on";
		my $title = "${prefix}Top-$top";
		$title =~ s/L/L\// if $title =~ /L[0-9]/;
		plot_contact_map("temp.rr", "${top}_$param_cont_type.png", $title, $file_log);
	}
	append2log("Show the maps");
	print "\n<table class=\"results\">";
	print "\n<th colspan=2 class=\"resultsheadingrow\"><b>Section ".Roman($sectionID++)."</b> Contact Maps (with native in grey background)</th>";
	print "\n<tr><td>";
	my $size = $html_display_width / 2 - 10;
	my $i = 1;
	foreach my $top (sort {$top_count{$a} <=> $top_count{$b}} keys %top_count){
		next if $top == "5";
		print "\n<a href=\"$web_root/cmap/${top}_$param_cont_type.png\" target=\"_blank\"><img src=\"$web_root/cmap/${top}_$param_cont_type.png\" width=\"$size\" align=\"middle\"/></a>";
		if($i == 1){
			print "\n</td><td>";
		}
		if($i == 3){
			print "\n</br>";
		}
		$i++;
		if($i > 1){
			$size = $html_display_width / 4 - 10;
		}
	}
	print "\n</td></tr>";
	print "\n</table>";
	print "\n</br>";
	chdir $dir_job or confess $!;
}

sub print_native_contact_maps{
	append2log("Plot native contact maps");
	my %selections = qw/all ALL/;
	foreach my $select (sort keys %selections){
		system_cmd("rm -f temp.rr");
		print2file("temp.rr", seq_rr($input_rr));
		my %rr = rr_rows_ordered_in_hash($native_rr, 100000, "all", $select);
		foreach my $row (sort {$a <=> $b} keys %rr){
			print2file("temp.rr", $rr{$row});
		}
		plot_cmap_native("temp.rr", "cmap_native_$select.png", "", $file_log);
	}
	append2log("Show the native contact maps");
	print "\n<table class=\"results\">";
	print "\n<th colspan=2 class=\"resultsheadingrow\"><b>Section ".Roman($sectionID++)."</b> Native's Contact Map</th>";
	print "\n<tr>";
	my $size = $html_display_width / 2 - 10;
	foreach my $select (sort keys %selections){
		print "\n<td><a href=\"$web_root/cmap_native_$select.png\" target=\"_blank\"><img src=\"$web_root/cmap_native_$select.png\" width=\"$size\" align=\"middle\"/></a></td>";
	}
	print "\n</tr>";
	print "\n</table>";
	print "\n</br>";
}

sub print_coord_num_rows{
	append2log("Calculating coordination number");
	foreach my $top (sort {$top_count{$a} <=> $top_count{$b}} keys %top_count){
		foreach my $src (sort keys %sources){
			my $count = $top_count{$top};
			$count = 100000 if $src eq "native";
			my $seq = undef;
			$seq = seq_rr($input_rr);
			# In case, we have the native pdb file, we want the reference sequence, to be that of pdb
			# So that R computations take that L into consideration
			$seq = seq_chain($native_pdb) if -f $native_pdb;
			my %rr = rr_rows_ordered_in_hash($input_rr, $count, $src, $param_cont_type);
			system_cmd("rm -f temp.rr");
			print2file("temp.rr", $seq);
			foreach my $row (sort {$a <=> $b} keys %rr){
				print2file("temp.rr", $rr{$row});
			}
			append2log("calculate coord_number for $src $top");
			$coord_number{$src}{$top} = rr_1D_coverage_line("temp.rr", $count);
		}
	}
	my $seq_num_row1 = "";
	my $seq_num_row2 = "";
	my $L = length(seq_rr($input_rr));
	foreach (my $i = 0; $i <= $L-3; $i = $i + 10){
		if ($i == 0){
			$seq_num_row1 .= sprintf "%-9s", "";
			$seq_num_row2 .= sprintf "%-9s", "";
		}
		else{
			$seq_num_row1 .= sprintf "%-10s", $i;
			$seq_num_row2 .= sprintf "%-10s", "|";
		}
	}
	append2log("Printing coordination number");
	print "\n<table class=\"results\">";
	print "\n<th colspan=13 class=\"resultsheadingrow\"><b>Section ".Roman($sectionID++)."</b> Coordination Numbers <font class=\"showdescription\"><a onclick=\"toggle_it('coornumdesc')\">see description</a></font></th>";
	# Hidden description
	print "\n<tr style=\"display:none\" id=\"coornumdesc\"><td colspan=100 class=\"resultdata\"><ul>";
	print "\n<li>A simple technique for 1D representation of predicted residue contacts is to assign numbers to each residue such that the numbers represent the number of contacts that the residue is involved in, also known as the coordination number.</li>";
	print "\n<li>For a 1D visualization by showing a single character number below the sequence, residues that are involved in less than 9 contacts are assigned numbers from 1 to 9, and the residues that are involved in more than 9 contacts are marked as '<b>*</b>'.</li>"; 
	print "\n<li>This visualization technique can show if contacts are clustered in a specific region or spread around evenly, and is effective when we have fewer contacts to analyze, for example L/10, L/5, L/2, L or even 2L. In addition, it is also convenient to compare contacts predicted by multiple sources.</li>";
	print "\n</ul></td></tr>";
	print "\n<tr>";
	# Results
	foreach my $top (sort {$top_count{$a} <=> $top_count{$b}} keys %top_count){
		my $prefix = "";
		$prefix = "Random-" if $param_na_rr_fl eq "on";
		my $id = "${prefix}Top-$top";
		$id =~ s/L/L\// if $id =~ /L[0-9]/;
		my $display = "";
		$display = "style=\"display:none\"" if $id ne "${prefix}Top-L/10";
		print "\n<tr><td><u><a onclick=\"toggle_it(\'$id\')\"><font color=\"blue\">".$id."</font></a></u></td></tr>";
		print "\n<tr $display id=\"$id\"><td class=\"coordinationnumber\" nowrap>";
		my $row_width = 100;
		my $num_of_breaks = 1 + int(length($rr_sequence)/$row_width);
		for(my $break = 0; $break < $num_of_breaks; $break++){
			my $start = $row_width * $break + 1;
			my $plus  = $row_width;
			if ($start + $plus > length($rr_sequence)){
				$plus = length($rr_sequence) - ($num_of_breaks-1) * $row_width;
			}
			my $this_seq_num_row1 = "   ".substr($seq_num_row1, $start + 2, $plus + 4);
			my $this_seq_num_row2 = "   ".substr($seq_num_row2, $start + 2, $plus + 4);
			if ($start + $plus > length($rr_sequence)){
				$this_seq_num_row1 = "   ".substr($seq_num_row1, $start + 2, $plus);
				$this_seq_num_row2 = "   ".substr($seq_num_row2, $start + 2, $plus);
			}
			$this_seq_num_row1 =~ s/\s/\&nbsp;/g;
			$this_seq_num_row2 =~ s/\s/\&nbsp;/g;
			print "$this_seq_num_row1</br>";
			print "$this_seq_num_row2</br>";
			for(my $i = $start - 1; $i < $start + $plus - 1; $i++){
				my $char = substr $rr_sequence, $i, 1;
				my $color = aa_color($char);
				printf "<font color=\"$color\">".$char."</font>", 0, 1;
			}
			print "\n<font color=\"gray\"> [sequence]</font></br>";
			if (defined $param_sec){
				for(my $i = $start - 1; $i < $start + $plus - 1; $i++){
					my $char = substr $param_sec, $i, 1;
					printf "<font color=\"red\">".$char."</font>", 0, 1   if $char eq "H";
					printf "<font color=\"blue\">".$char."</font>", 0, 1  if $char eq "E";
					printf "<font color=\"green\">".$char."</font>", 0 ,1 if $char eq "C";
				}
				print "\n<font color=\"gray\"> [secondary structure]</font></br>";
			}
			foreach my $src (sort keys %sources){
				my $coverage_line = $coord_number{$src}{$top};
				my @C = split /\s\[/, $coverage_line;
				print "".substr($C[0], $start - 1, $plus)."<font color=\"gray\"> [$src: ".$C[1]."</font>";
				print "\n</br>";
			}
			print "\n</br>";
		}
		print "\n</td></tr>";
	}
	print "\n</table>";
	print "\n</br>";
}

sub rrhash2chimerascript{
	my $rrhash = shift;
	my $chimera_script = shift;
	my %rr = %{$rrhash};
	confess "Empty native atoms file" if not scalar keys %min_true_d_atoms;
	print2file($chimera_script, "# How to visualize these contacts in UCSF Chimera?");
	print2file($chimera_script, "# 1. Download and install UCSF Chimera 1.10 or later if you do not have it already.");
	print2file($chimera_script, "# 2. Save this file, say, 'all_sr.com'.");
	print2file($chimera_script, "# 3. Use 'File->Open' in UCSF Chimera to open 'all_sr.com'.");
	print2file($chimera_script, "");
	print2file($chimera_script, "background solid white;");
	print2file($chimera_script, "open $web_root/native.pdb;");
	print2file($chimera_script, "focus;");
	print2file($chimera_script, "rainbow;");
	print2file($chimera_script, "labelopt resinfo \%(number)s");
	print2file($chimera_script, "represent wire;");
	print2file($chimera_script, "show;");
	print2file($chimera_script, "rlabel;");
	print2file($chimera_script, "focus;");
	foreach (sort {$a <=> $b} keys %rr){
		my @C = split /\s+/, $rr{$_};
		my $atoms = $min_true_d_atoms{$C[0]." ".$C[1]};
		my @A = split /\s+/, $atoms;
		confess "Atoms must be defined for rr row $_ [".$C[0]." ".$C[1]."]" if not $A[1];
		print2file($chimera_script, "distance :".$C[0]."@".$A[0]." :".$C[1]."@".$A[1]);
	}
	print2file($chimera_script, "setattr g lineWidth 2;");
	print2file($chimera_script, "setattr g color black;");
	print2file($chimera_script, "turn z 2 360 precess 20;");
}

sub print_precision{
	append2log("Calculate precision");
	foreach my $top (sort {$top_count_order{$a} <=> $top_count_order{$b}} keys %top_count_order){
		foreach my $src (sort keys %sources){
			my %rr_selected = %{$selected_contacts{$src." ".$top}};
			confess "\nERROR! Empty RR file! Not expecting!" if not scalar keys %rr_selected;
			$precision{$src}{$top} = calc_precision(\%rr_selected);
		}
	}
	append2log("Plot precision graphs");
	foreach my $s(sort keys %sources){
		foreach my $t (sort {$top_count_order{$a} <=> $top_count_order{$b}} keys %top_count_order){
			print2file("precision_$param_cont_type.txt", "$s\t".$top_count_order{$t}."\t".$precision{$s}{$t}."\t$t");
		}
	}
	print2file("precision_$param_cont_type.R", "png(\"precision_$param_cont_type.png\", width=2500, height=1400, res=300)");
	print2file("precision_$param_cont_type.R", "par(mar=c(4.5,4.5,0.2,11), xpd = TRUE)");
	print2file("precision_$param_cont_type.R", "infile = \"precision_$param_cont_type.txt\"");
	print2file("precision_$param_cont_type.R", "original <- read.table(file=infile, header=F)");
	print2file("precision_$param_cont_type.R", "sorted <- original[order(original\$V1, original\$V2),]");
	print2file("precision_$param_cont_type.R", "y_range <- range(0, original\$V3)");
	print2file("precision_$param_cont_type.R", "sources <- unique(original\$V1)");
	print2file("precision_$param_cont_type.R", "start_flag = 1");
	print2file("precision_$param_cont_type.R", "src_id = 1");
	print2file("precision_$param_cont_type.R", "for(src in sources){");
	print2file("precision_$param_cont_type.R", "\tselected <- original\$V3[which(original\$V1 == src)]");
	print2file("precision_$param_cont_type.R", "\tif(start_flag == 1){");
	print2file("precision_$param_cont_type.R", "\t\tplot(selected, col=src_id, type=\"l\", lwd = 2, xaxt=\'n\', ylim=y_range, ylab=\"Precision\", xlab=\"Top Ranked Contacts\")");
	print2file("precision_$param_cont_type.R", "\t\tstart_flag = 0");
	print2file("precision_$param_cont_type.R", "\t}");
	print2file("precision_$param_cont_type.R", "\telse{");
	print2file("precision_$param_cont_type.R", "\t\tlines(selected, col=src_id, lty=src_id,  lwd = 3)");
	print2file("precision_$param_cont_type.R", "\t}");
	print2file("precision_$param_cont_type.R", "\tsrc_id = src_id + 1");
	print2file("precision_$param_cont_type.R", "}");
	print2file("precision_$param_cont_type.R", "axis(1, at=1:6, labels=c(\"Top5\",\"L/10\",\"L/5\",\"L/2\",\"L\",\"2L\"))");
	print2file("precision_$param_cont_type.R", "if (length(sources) > 1){");
	print2file("precision_$param_cont_type.R", "\tlegend(\"topright\", inset=c(-0.4,0), legend = sources, cex=1.0, lwd = 3, col = seq(1, length(sources), by = 1), lty = seq(1, length(sources), by = 1))");
	print2file("precision_$param_cont_type.R", "}");
	system_cmd("Rscript precision_$param_cont_type.R");
	append2log("Show precision");

	print "\n<table class=\"results\">";
	print "\n<th colspan=13 class=\"resultsheadingrow\"><b>Section ".Roman($sectionID++)."</b> Precision <font class=\"showdescription\"><a onclick=\"toggle_it('precisiondesc')\">see description</a></font></th>";
	# Hidden description
	print "\n<tr style=\"display:none\" id=\"precisiondesc\"><td colspan=100 class=\"resultdata\"><ul>";
	print "\n<li>Precision = TP/(TP+FP), where TP = True Positives, FP = False Positives, is calculated as the ratio of the number of correctly predicted contacts and the total number of predicted contacts.</li>";
	print "\n</ul></td></tr>";
	print "\n<tr>";
	# Results
	print "\n<tr><td>";
	# a sub table of precision
	print "\n<table class=\"smalldatabox\">";
	print "\n<th class=\"resultsdataheadingrow\">RR-File</th>";
	foreach my $top (sort {$top_count_order{$a} <=> $top_count_order{$b}} keys %top_count_order){
		$top =~ s/L/L\// if $top =~ /L[0-9]/;
		print "\n<th class=\"resultsdataheadingrow\">Top $top</th>";
	}
	foreach my $src (sort keys %sources){
		print "\n<tr>";
		print "\n<td class=\"smalldataboxvalues\">$src</td>";
		foreach my $top (sort {$top_count_order{$a} <=> $top_count_order{$b}} keys %top_count_order){
			print "\n<td class=\"smalldataboxvalues\">".$precision{$src}{$top}."</td>";
		}
		print "\n</tr>";
	}
	# close sub table
	print "\n</table>";
	print "\n</td><td>";
	print "\n<a href=\"$web_root/precision_$param_cont_type.png\" target=\"_blank\"><img src=\"$web_root/precision_$param_cont_type.png\" width=\"".($html_display_width / 2 - 100)."\" align=\"middle\"/></a>";
	print "\n</td></tr></table>";
	print "\n</br>";
}

sub print_coverage{
	append2log("Calculate coverage");
	foreach my $top (keys %top_count_order){
		foreach my $src (sort keys %sources){
			my %rr_selected = %{$selected_contacts{$src." ".$top}};
			confess "\nERROR! Empty RR file! Not expecting this!" if not scalar keys %rr_selected;
			$coverage{$src}{$top} = calc_coverage(\%rr_selected);
		}
	}
	append2log("Plot coverage");
	foreach my $s(sort keys %sources){
		foreach my $t (sort {$top_count_order{$a} <=> $top_count_order{$b}} keys %top_count_order){
			print2file("coverage_$param_cont_type.txt", "$s\t".$top_count_order{$t}."\t".$coverage{$s}{$t}."\t$t");
		}
	}
	system_cmd("rm -f coverage_$param_cont_type.R");
	print2file("coverage_$param_cont_type.R", "png(\"coverage_$param_cont_type.png\", width=2500, height=1400, res=300)");
	print2file("coverage_$param_cont_type.R", "par(mar=c(4.5,4.5,0.2,11), xpd = TRUE)");
	print2file("coverage_$param_cont_type.R", "infile = \"coverage_$param_cont_type.txt\"");
	print2file("coverage_$param_cont_type.R", "original <- read.table(file=infile, header=F)");
	print2file("coverage_$param_cont_type.R", "sorted <- original[order(original\$V1, original\$V2),]");
	print2file("coverage_$param_cont_type.R", "y_range <- range(0, original\$V3)");
	print2file("coverage_$param_cont_type.R", "sources <- unique(original\$V1)");
	print2file("coverage_$param_cont_type.R", "start_flag = 1");
	print2file("coverage_$param_cont_type.R", "src_id = 1");
	print2file("coverage_$param_cont_type.R", "for(src in sources){");
	print2file("coverage_$param_cont_type.R", "\tselected <- original\$V3[which(original\$V1 == src)]");
	print2file("coverage_$param_cont_type.R", "\tif(start_flag == 1){");
	print2file("coverage_$param_cont_type.R", "\t\tplot(selected, col=src_id, type=\"l\", lwd = 2, xaxt=\'n\', ylim=y_range, ylab=\"Coverage\", xlab=\"Top Ranked Contacts\")");
	print2file("coverage_$param_cont_type.R", "\t\tstart_flag = 0");
	print2file("coverage_$param_cont_type.R", "\t}");
	print2file("coverage_$param_cont_type.R", "\telse{");
	print2file("coverage_$param_cont_type.R", "\t\tlines(selected, col=src_id, lty=src_id,  lwd = 3)");
	print2file("coverage_$param_cont_type.R", "\t}");
	print2file("coverage_$param_cont_type.R", "\tsrc_id = src_id + 1");
	print2file("coverage_$param_cont_type.R", "}");
	print2file("coverage_$param_cont_type.R", "axis(1, at=1:6, labels=c(\"Top5\",\"L/10\",\"L/5\",\"L/2\",\"L\",\"2L\"))");
	print2file("coverage_$param_cont_type.R", "if (length(sources) > 1){");
	print2file("coverage_$param_cont_type.R", "\tlegend(\"topright\", inset=c(-0.4,0), legend = sources, cex=1.0, lwd = 3, col = seq(1, length(sources), by = 1), lty = seq(1, length(sources), by = 1))");
	print2file("coverage_$param_cont_type.R", "}");
	system_cmd("Rscript coverage_$param_cont_type.R");
	append2log("Show coverage");
	print "\n<table class=\"results\">";
	print "\n<th colspan=13 class=\"resultsheadingrow\"><b>Section ".Roman($sectionID++)."</b> Coverage <font class=\"showdescription\"><a onclick=\"toggle_it('coveragedesc')\">see description</a></font></th>";
	# Hidden description
	print "\n<tr style=\"display:none\" id=\"coveragedesc\"><td colspan=100 class=\"resultdata\"><ul>";
	print "\n<li>Coverage = 100 * TP/N, where, TP is 'Number of true positive contacts' and N is 'Number of native contacts'. Same contact definition is used for calculating both numerator and denominator.</li>";
	print "\n</ul></td></tr>";
	print "\n<tr>";
	# Results
	print "\n<tr><td>";
	# a sub table for precise data
	print "\n<table class=\"smalldatabox\">";
	print "\n<th class=\"resultsdataheadingrow\">RR-File</td>";
	foreach my $top (sort {$top_count_order{$a} <=> $top_count_order{$b}} keys %top_count_order){
		$top =~ s/L/L\// if $top =~ /L[0-9]/;
		my $prefix = "";
		$prefix = "Random" if $param_na_rr_fl;
		print "\n<th class=\"resultsdataheadingrow\">$prefix Top $top</th>";
	}
	foreach my $src (sort keys %sources){
		print "\n<tr><td class=\"smalldataboxvalues\">$src</td>";
		foreach my $top (sort {$top_count_order{$a} <=> $top_count_order{$b}} keys %top_count_order){
			print "\n<td class=\"smalldataboxvalues\">".$coverage{$src}{$top}."</td>";
		}
		print "\n</tr>";
	}
	# close sub table
	print "\n</table>"; 
	print "\n</td><td>";
	print "\n<a href=\"$web_root/coverage_$param_cont_type.png\" target=\"_blank\"><img src=\"$web_root/coverage_$param_cont_type.png\" width=\"".($html_display_width / 2 - 100)."\" align=\"middle\"/></a>";
	print "\n</td></tr></table>";
	print "\n</br>";
}

sub print_xd{
	append2log("Calculate Xd");
	my $seq = undef;
	foreach my $top (sort {$top_count_order{$a} <=> $top_count_order{$b}} keys %top_count_order){
		foreach my $src (sort keys %sources){
			$seq = seq_rr($input_rr);
			# In case, we have the native pdb file, we want the reference sequence, to be that of pdb
			# So that R computations take that L into consideration
			$seq = seq_chain($native_pdb) if -f $native_pdb;
			append2log("RR Select: $top");
			append2log("RR source: $src");
			my %rr_selected = rr_rows_ordered_in_hash($input_rr, $top_count{$top}, $src, $param_cont_type);
			$xd{$src}{$top} = calc_xd(\%rr_selected);
		}
	}
	append2log("Plot Xd graphs");
	foreach my $s(sort keys %sources){
		foreach my $t (sort {$top_count_order{$a} <=> $top_count_order{$b}} keys %top_count_order){
			print2file("xd_$param_cont_type.txt", "$s\t".$top_count_order{$t}."\t".$xd{$s}{$t}."\t$t");
		}
	}
	print2file("xd_$param_cont_type.R", "png(\"xd_$param_cont_type.png\", width=2500, height=1400, res=300)");
	print2file("xd_$param_cont_type.R", "par(mar=c(4.5,4.5,0.3,11), xpd = TRUE)");
	print2file("xd_$param_cont_type.R", "infile = \"xd_$param_cont_type.txt\"");
	print2file("xd_$param_cont_type.R", "original <- read.table(file=infile, header=F)");
	print2file("xd_$param_cont_type.R", "sorted <- original[order(original\$V1, original\$V2),]");
	print2file("xd_$param_cont_type.R", "y_range <- range(0, original\$V3)");
	print2file("xd_$param_cont_type.R", "sources <- unique(original\$V1)");
	print2file("xd_$param_cont_type.R", "start_flag = 1");
	print2file("xd_$param_cont_type.R", "src_id = 1");
	print2file("xd_$param_cont_type.R", "for(src in sources){");
	print2file("xd_$param_cont_type.R", "\tselected <- original\$V3[which(original\$V1 == src)]");
	print2file("xd_$param_cont_type.R", "\tif(start_flag == 1){");
	print2file("xd_$param_cont_type.R", "\t\tplot(selected, col=src_id, type=\"l\", lwd = 2, xaxt=\'n\', ylim=y_range, ylab=\"Xd\", xlab=\"Top Ranked Contacts\")");
	print2file("xd_$param_cont_type.R", "\t\tstart_flag = 0");
	print2file("xd_$param_cont_type.R", "\t}");
	print2file("xd_$param_cont_type.R", "\telse{");
	print2file("xd_$param_cont_type.R", "\t\tlines(selected, col=src_id, lty=src_id,  lwd = 3)");
	print2file("xd_$param_cont_type.R", "\t}");
	print2file("xd_$param_cont_type.R", "\tsrc_id = src_id + 1");
	print2file("xd_$param_cont_type.R", "}");
	print2file("xd_$param_cont_type.R", "axis(1, at=1:6, labels=c(\"Top5\",\"L/10\",\"L/5\",\"L/2\",\"L\",\"2L\"))");
	print2file("xd_$param_cont_type.R", "if (length(sources) > 1){");
	print2file("xd_$param_cont_type.R", "\tlegend(\"topright\", inset=c(-0.4,0), legend = sources, cex=1.0, lwd = 3, col = seq(1, length(sources), by = 1), lty = seq(1, length(sources), by = 1))");
	print2file("xd_$param_cont_type.R", "}");
	system_cmd("Rscript xd_$param_cont_type.R");
	append2log("Show Xd"); 
	print "\n<table class=\"results\">";
	print "\n<th colspan=13 class=\"resultsheadingrow\"><b>Section ".Roman($sectionID++)."</b> X<sub>d</sub> (Distance distribution) <font class=\"showdescription\"><a onclick=\"toggle_it('xddesc')\">see description</a></font></th>";
	# Hidden description
	print "\n<tr style=\"display:none\" id=\"xddesc\"><td colspan=100 class=\"resultdata\"><ul>";
	print "\n<li>X<sub>d</sub> or Distance distribution = (P<sub>p1</sub>-P<sub>a1</sub>)/(15 * d<sib>1</sub>) + .. + (P<sub>pi</sub>-P<sub>ai</sub>)/(15 * d<sub>i</sub>)+ .. + (P<sub>p15</sub>-P<sub>a15</sub>)/(15 * d<sib>15</sub>), where P<sub>pi</sub> is the fraction of predicted contacts in distance bin i, P<sub>ai</sub> is the fraction of all residue pairs in distance bin i.</li>";
	print "\n<li>Xd measures the measures the weighted harmonic average difference between the distance distribution of predicted contacts and the all-pairs.</li>"; 
	print "\n<li>There are 15 distance bins covering the range from 0 to 60 Å. The 15 bins include ranges of distances from 0 to 4 Å, 4 Å to 8 Å, 8 Å to 12 Å, etc.</li>";
	print "\n<li>This score estimates the deviation of the distribution of distances in the list of contacts from the distribution of distances in all pairs of residues in the protein. Since CASP6 precision and Xd have been consistently used for contact evaluation in all the CASP competitions afterwards, including the most recent CASP10 competition.</li>";
	print "\n</ul></td></tr>";
	print "\n<tr>";
	# Results
	print "\n<tr><td>";
	print "\n<a href=\"$web_root/xd_$param_cont_type.png\" target=\"_blank\"><img src=\"$web_root/xd_$param_cont_type.png\" width=\"".($html_display_width / 2 - 100)."\" align=\"middle\"/></a>";
	print "\n</td><td>";
	# a sub table of values
	print "\n<table class=\"smalldatabox\">";
	print "\n<th class=\"resultsdataheadingrow\">RR-File</td>";
	foreach my $top (sort {$top_count_order{$a} <=> $top_count_order{$b}} keys %top_count_order){
		$top =~ s/L/L\// if $top =~ /L[0-9]/;
		my $prefix = "";
		$prefix = "Random" if $param_na_rr_fl;
		print "\n<th class=\"resultsdataheadingrow\">$prefix Top $top</th>";
	}
	foreach my $src (sort keys %sources){
		print "\n<tr><td class=\"smalldataboxvalues\">$src</td>";
		foreach my $top (sort {$top_count_order{$a} <=> $top_count_order{$b}} keys %top_count_order){
			print "\n<td class=\"smalldataboxvalues\">".$xd{$src}{$top}."</td>";
		}
		print "\n</tr>";
	}
	# close sub table
	print "\n</table>";
	print "\n</td></tr></table>";
	print "\n</br>";
}

sub print_spread{
	append2log("Calculate Spread");
	my $seq = undef;
	foreach my $top (sort {$top_count_order{$a} <=> $top_count_order{$b}} keys %top_count_order){
		foreach my $src (sort keys %sources){
			$seq = seq_rr($input_rr);
			# In case, we have the native pdb file, we want the reference sequence, to be that of pdb
			# So that R computations take that L into consideration
			$seq = seq_chain($native_pdb) if -f $native_pdb;
			append2log("RR-Select: $top RR-source: $src");
			my %rr_selected = rr_rows_ordered_in_hash($input_rr, $top_count{$top}, $src, $param_cont_type);
			system_cmd("rm -f temp.rr");
			print2file("temp.rr", $seq);
			foreach (sort {$a <=> $b} keys %rr_selected){
				print2file("temp.rr", $rr_selected{$_});
			}
			$spread{$src}{$top} = calc_spread("temp.rr");
		}
	}
	append2log("Plot Spread graphs");
	foreach my $s(sort keys %sources){
		foreach my $t (sort {$top_count_order{$a} <=> $top_count_order{$b}} keys %top_count_order){
			print2file("spread_$param_cont_type.txt", "$s\t".$top_count_order{$t}."\t".$spread{$s}{$t}."\t$t");
		}
	}
	print2file("spread_$param_cont_type.R", "png(\"spread_$param_cont_type.png\", width=2500, height=1400, res=300)");
	print2file("spread_$param_cont_type.R", "par(mar=c(4.5,4.5,0.2,11), xpd = TRUE)");
	print2file("spread_$param_cont_type.R", "infile = \"spread_$param_cont_type.txt\"");
	print2file("spread_$param_cont_type.R", "original <- read.table(file=infile, header=F)");
	print2file("spread_$param_cont_type.R", "sorted <- original[order(original\$V1, original\$V2),]");
	print2file("spread_$param_cont_type.R", "y_range <- range(0, original\$V3)");
	print2file("spread_$param_cont_type.R", "sources <- unique(original\$V1)");
	print2file("spread_$param_cont_type.R", "start_flag = 1");
	print2file("spread_$param_cont_type.R", "src_id = 1");
	print2file("spread_$param_cont_type.R", "for(src in sources){");
	print2file("spread_$param_cont_type.R", "\tselected <- original\$V3[which(original\$V1 == src)]");
	print2file("spread_$param_cont_type.R", "\tif(start_flag == 1){");
	print2file("spread_$param_cont_type.R", "\t\tplot(selected, col=src_id, type=\"l\", lwd = 2, xaxt=\'n\', ylim=y_range, ylab=\"Spread\", xlab=\"Top Ranked Contacts\")");
	print2file("spread_$param_cont_type.R", "\t\tstart_flag = 0");
	print2file("spread_$param_cont_type.R", "\t}");
	print2file("spread_$param_cont_type.R", "\telse{");
	print2file("spread_$param_cont_type.R", "\t\tlines(selected, col=src_id, lty=src_id,  lwd = 3)");
	print2file("spread_$param_cont_type.R", "\t}");
	print2file("spread_$param_cont_type.R", "\tsrc_id = src_id + 1");
	print2file("spread_$param_cont_type.R", "}");
	print2file("spread_$param_cont_type.R", "axis(1, at=1:6, labels=c(\"Top5\",\"L/10\",\"L/5\",\"L/2\",\"L\",\"2L\"))");
	print2file("spread_$param_cont_type.R", "if (length(sources) > 1){");
	print2file("spread_$param_cont_type.R", "\tlegend(\"topright\", inset=c(-0.4,0), legend = sources, cex=1.0, lwd = 3, col = seq(1, length(sources), by = 1), lty = seq(1, length(sources), by = 1))");
	print2file("spread_$param_cont_type.R", "}");
	system_cmd("Rscript spread_$param_cont_type.R");
	append2log("Show Spread"); 
	
	print "\n<table class=\"results\">";
	print "\n<th colspan=13 class=\"resultsheadingrow\"><b>Section ".Roman($sectionID++)."</b> Spread <font class=\"showdescription\"><a onclick=\"toggle_it('spreaddesc')\">see description</a></font></th>";
	# Hidden description
	print "\n<tr style=\"display:none\" id=\"spreaddesc\"><td colspan=100 class=\"resultdata\"><ul>";
	print "\n<li>Spread is computed using contact maps. It is computed as the mean of the distances from every experimental (crystal structure) contact to the nearest predicted contact in 2D contact map.</li>";
	print "\n</ul></td></tr>";
	print "\n<tr>";
	# Results
	print "\n<tr><td>";
	# a sub table for precise data
	print "\n<table class=\"smalldatabox\">";
	print "\n<th class=\"resultsdataheadingrow\">RR-File</td>";
	foreach my $top (sort {$top_count_order{$a} <=> $top_count_order{$b}} keys %top_count_order){
		$top =~ s/L/L\// if $top =~ /L[0-9]/;
		my $prefix = "";
		$prefix = "Random" if $param_na_rr_fl;
		print "\n<th class=\"resultsdataheadingrow\">$prefix Top $top</th>";
	}
	print "\n</tr>";
	foreach my $src (sort keys %sources){
		print "\n<tr><td class=\"smalldataboxvalues\">$src</td>";
		foreach my $top (sort {$top_count_order{$a} <=> $top_count_order{$b}} keys %top_count_order){
			print "\n<td class=\"smalldataboxvalues\">".$spread{$src}{$top}."</td>";
		}
		print "\n</tr>";
	}
	# close sub table
	print "\n</table>"; 
	print "\n</td><td>";
	print "\n<a href=\"$web_root/spread_$param_cont_type.png\" target=\"_blank\"><img src=\"$web_root/spread_$param_cont_type.png\" width=\"".($html_display_width / 2 - 100)."\" align=\"middle\"/></a>";
	print "\n</td></tr></table>";
	print "\n</br>";
}

sub print_FP_assessment{
	append2log("Calculate False Positive Error");
	my $seq = undef;
	foreach my $top (sort {$top_count_order{$a} <=> $top_count_order{$b}} keys %top_count_order){
		foreach my $src (sort keys %sources){
			my %rr_selected = %{$selected_contacts{$src." ".$top}};
			confess "??" if not scalar keys %rr_selected;
			$fp{$src}{$top} = calc_fp(\%rr_selected);
		}
	}
	append2log("Plot fp graphs");
	foreach my $s(sort keys %sources){
		foreach my $t (sort {$top_count_order{$a} <=> $top_count_order{$b}} keys %top_count_order){
			print2file("fp_$param_cont_type.txt", "$s\t".$top_count_order{$t}."\t".$fp{$s}{$t}."\t$t");
		}
	}
	print2file("fp_$param_cont_type.R", "png(\"fp_$param_cont_type.png\", width=2500, height=1400, res=300)");
	print2file("fp_$param_cont_type.R", "par(mar=c(4.5,4.5,0.2,11), xpd = TRUE)");
	print2file("fp_$param_cont_type.R", "infile = \"fp_$param_cont_type.txt\"");
	print2file("fp_$param_cont_type.R", "original <- read.table(file=infile, header=F)");
	print2file("fp_$param_cont_type.R", "sorted <- original[order(original\$V1, original\$V2),]");
	print2file("fp_$param_cont_type.R", "y_range <- range(0, original\$V3)");
	print2file("fp_$param_cont_type.R", "sources <- unique(original\$V1)");
	print2file("fp_$param_cont_type.R", "start_flag = 1");
	print2file("fp_$param_cont_type.R", "src_id = 1");
	print2file("fp_$param_cont_type.R", "for(src in sources){");
	print2file("fp_$param_cont_type.R", "\tselected <- original\$V3[which(original\$V1 == src)]");
	print2file("fp_$param_cont_type.R", "\tif(start_flag == 1){");
	print2file("fp_$param_cont_type.R", "\t\tplot(selected, col=src_id, type=\"l\", lwd = 2, xaxt=\'n\', ylim=y_range, ylab=\"Mean FP Error\", xlab=\"Top Ranked Contacts\")");
	print2file("fp_$param_cont_type.R", "\t\tstart_flag = 0");
	print2file("fp_$param_cont_type.R", "\t}");
	print2file("fp_$param_cont_type.R", "\telse{");
	print2file("fp_$param_cont_type.R", "\t\tlines(selected, col=src_id, lty=src_id,  lwd = 3)");
	print2file("fp_$param_cont_type.R", "\t}");
	print2file("fp_$param_cont_type.R", "\tsrc_id = src_id + 1");
	print2file("fp_$param_cont_type.R", "}");
	print2file("fp_$param_cont_type.R", "axis(1, at=1:6, labels=c(\"Top5\",\"L/10\",\"L/5\",\"L/2\",\"L\",\"2L\"))");
	print2file("fp_$param_cont_type.R", "if (length(sources) > 1){");
	print2file("fp_$param_cont_type.R", "\tlegend(\"bottomright\", inset=c(-0.4,0), legend = sources, cex=1.0, lwd = 3, col = seq(1, length(sources), by = 1), lty = seq(1, length(sources), by = 1))");
	print2file("fp_$param_cont_type.R", "}");
	system_cmd("Rscript fp_$param_cont_type.R");
	append2log("Show Xd"); 
	print "\n<table class=\"results\">";
	print "\n<th colspan=13 class=\"resultsheadingrow\"><b>Section ".Roman($sectionID++)."</b> Mean False Positive Error <font class=\"showdescription\"><a onclick=\"toggle_it('mfpedesc')\">see description</a></font></th>";
	# Hidden description
	print "\n<tr style=\"display:none\" id=\"mfpedesc\"><td colspan=100 class=\"resultdata\"><ul>";
	print "\n<li>Mean FP Error = avg(d<sub>error</sub>); d<sub>error</sub> = 0 if the predicted contact is correct. Otherwise, d<sub>error</sub> = d - T, where d is the actual distance of a false positive contact in native structure and T is the distance threshold for the predicted contact.</li>";
	print "\n</ul></td></tr>";
	print "\n<tr>";
	# Results
	print "\n<tr><td>";
	print "\n<a href=\"$web_root/fp_$param_cont_type.png\" target=\"_blank\"><img src=\"$web_root/fp_$param_cont_type.png\" width=\"".($html_display_width / 2 - 100)."\" align=\"middle\"/></a>";
	print "\n</td><td>";
	# a sub table of values
	print "\n<table class=\"smalldatabox\">";
	print "\n<th class=\"resultsdataheadingrow\">RR-File</td>";
	foreach my $top (sort {$top_count_order{$a} <=> $top_count_order{$b}} keys %top_count_order){
		$top =~ s/L/L\// if $top =~ /L[0-9]/;
		print "\n<th class=\"resultsdataheadingrow\">Top $top</th>";
	}
	foreach my $src (sort keys %sources){
		print "\n<tr><td class=\"smalldataboxvalues\">$src</td>";
		foreach my $top (sort {$top_count_order{$a} <=> $top_count_order{$b}} keys %top_count_order){
			print "\n<td class=\"smalldataboxvalues\">".$fp{$src}{$top}."</td>";
		}
		print "\n</tr>";
	}
	# close sub table
	print "\n</table>";
	print "\n</td></tr></table>";
	print "\n</br>";
}

sub dssp_result{
	my $file_pdb = shift;
	my $selection = shift;
	confess "ERROR! $file_pdb does not exist!" if not -f $file_pdb;
	confess "ERROR! selection not defined!" if not defined $selection;
	my %RESIDUE = ();
	my %SS  = ();
	my %PHI = ();
	my %PSI = ();
	my $command = "$program_dssp $file_pdb | grep -C 0 -A 1000 \'  #  RESIDUE\' | tail -n +2";
	my @dssp_rows = `$command`;
	foreach(@dssp_rows){
		my $rnum = substr($_,  5, 6);
		my $res  = substr($_, 13, 1);
		my $sstr = substr($_, 16, 1);
		my $phia = substr($_,103, 6);
		my $psia = substr($_,109, 6);
		$rnum =~ s/\s+//g;
		$res  =~ s/\s+//g;
		$sstr =~ s/\s+//g;
		$phia =~ s/\s+//g;
		$psia =~ s/\s+//g;
		# alternate residue locations may have alphabets
		#$rnum =~ s/[^0-9]//g;
		$res  =~ s/[^A-Z]//g;
		$sstr =~ s/[^A-Z]//g;
		next if length($rnum) < 1;
		confess ":( residue not defined for $rnum" if length($res) < 1;
		confess ":( phi not defined for $rnum" if length($phia) < 1;
		confess ":( psi not defined for $rnum" if length($psia) < 1;
		$sstr = "C" if length($sstr) < 1;
		$sstr =~ s/\./C/g;
		$sstr =~ s/S/C/g;
		$sstr =~ s/T/C/g;
		$sstr =~ s/B/C/g;
		$sstr =~ s/G/C/g;
		$RESIDUE{$rnum} = $res;
		$SS{$rnum}      = $sstr;
		$PHI{$rnum}     = $phia;
		$PSI{$rnum}     = $psia;
	}

	# special temporary change
	return if not scalar keys %PHI;

	confess "$file_pdb seems empty!" if not scalar keys %PHI;
	return %SS if ($selection eq "ss");
	return %PHI if ($selection eq "phi");
	return %PSI if ($selection eq "psi");
	confess "ERROR! Invalid selection string!";
}

sub dssp_ss_pdb_to_file{
	my $file_pdb = shift;
	my $file_ss = shift;
	confess ":(" if not -f $file_pdb;
	my %ss = dssp_result($file_pdb, "ss");
	foreach(sort {$a <=> $b} keys %ss){
		print2file($file_ss, $_." ".$ss{$_});
	}
	confess "looks like DSSP failed!" if not -s $file_ss;
}

sub plot_contact_map{
	my $file_rr = shift;
	my $output_png = shift;
	my $id = shift;
	append2log("Plotting $id");
	system_cmd("rm -f cmap.R");
	print2file("cmap.R", "png(\"$output_png\", width=1600, height=1600, res=300)");
	print2file("cmap.R", "infile = \"$file_rr\"");
	print2file("cmap.R", "con <- file(infile, \"rt\")"); 
	print2file("cmap.R", "seqlen <- readLines(con, 1)");
	print2file("cmap.R", "L <- nchar(seqlen)");
	print2file("cmap.R", "data <- read.table(infile, skip = 1, header = FALSE)");
	print2file("cmap.R", "col_count <- ncol(data)");
	print2file("cmap.R", "data\$V7 <- abs(data\$V1-data\$V2)");
	print2file("cmap.R", "data\$V7[data\$V7 > 23] <- 200");
	print2file("cmap.R", "data\$V7[data\$V7 < 12] <- 100");
	print2file("cmap.R", "data\$V7[data\$V7 < 24] <- 150");
	print2file("cmap.R", "data\$V7[data\$V7 == 100] <- 4");
	print2file("cmap.R", "data\$V7[data\$V7 == 150] <- 3");
	print2file("cmap.R", "data\$V7[data\$V7 == 200] <- 2");
	print2file("cmap.R", "data\$V8[data\$V7 == 4] <- \"Short-Range\"");
	print2file("cmap.R", "data\$V8[data\$V7 == 3] <- \"Medium-Range\"");
	print2file("cmap.R", "data\$V8[data\$V7 == 2] <- \"Long-Range\"");
	print2file("cmap.R", "data\$V9 = as.numeric(data\$V6)");
	print2file("cmap.R", "par(mar=c(2,2,0.5,0.5))");
	print2file("cmap.R", "# Plot for an RR file with multiple sources");
	print2file("cmap.R", "native <- subset(data, data\$V6 == \"PDB\")");
	print2file("cmap.R", "plot(native\$V1, native\$V2, col=rgb(0,0,0,0.1), pch=0, xlab = NULL, ylab = NULL, ylim=c(L, 1), xlim=c(1, L))");
	print2file("cmap.R", "predic <- subset(data, data\$V6 != \"PDB\")");
	print2file("cmap.R", "points(predic\$V1, predic\$V2, col=predic\$V6, xlab = NULL, ylab = NULL, ylim=c(L, 1), xlim=c(1, L), pch = predic\$V9)");
	print2file("cmap.R", "legend(\"topright\", bty=\"n\", legend=unique(predic\$V6), pch = unique(predic\$V9), col=unique(predic\$V9))");
	print2file("cmap.R", "text(0, 0, \"$id\", cex=2.0, adj = c(0,1))");
	#print2file("cmap.R", "dev.off()");
	system_cmd("Rscript cmap.R");
}

sub plot_cmap_native{
	my $file_rr = shift;
	my $output_png = shift;
	my $id = shift;

	# Calculate Native SS and print to file
	# Throws errors for cases like http://www.rcsb.org/pdb/files/1aap.pdb
	# dssp_ss_pdb_to_file($native_pdb, "native.ss");

	append2log("Plotting $id");
	system_cmd("rm -f cmap.R");
	print2file("cmap.R", "png(\"$output_png\", width=1600, height=1600, res=300)");
	print2file("cmap.R", "infile = \"$file_rr\"");
	print2file("cmap.R", "con <- file(infile, \"rt\")"); 
	print2file("cmap.R", "seqlen <- readLines(con, 1)");
	print2file("cmap.R", "L <- nchar(seqlen)");
	print2file("cmap.R", "data <- read.table(infile, skip = 1, header = FALSE)");
	print2file("cmap.R", "col_count <- ncol(data)");
	print2file("cmap.R", "data\$V7 <- abs(data\$V1-data\$V2)");
	print2file("cmap.R", "data\$V7[data\$V7 > 23] <- 200");
	print2file("cmap.R", "data\$V7[data\$V7 < 12] <- 100");
	print2file("cmap.R", "data\$V7[data\$V7 < 24] <- 150");
	print2file("cmap.R", "data\$V7[data\$V7 == 100] <- 4");
	print2file("cmap.R", "data\$V7[data\$V7 == 150] <- 3");
	print2file("cmap.R", "data\$V7[data\$V7 == 200] <- 2");
	print2file("cmap.R", "data\$V8[data\$V7 == 4] <- \"Short-Range\"");
	print2file("cmap.R", "data\$V8[data\$V7 == 3] <- \"Medium-Range\"");
	print2file("cmap.R", "data\$V8[data\$V7 == 2] <- \"Long-Range\"");
	print2file("cmap.R", "par(mar=c(2,2,0.5,0.5))");
	print2file("cmap.R", "if(col_count == 6){");
	print2file("cmap.R", "\t# Plot for an RR file with multiple sources");
	print2file("cmap.R", "\tplot(data\$V1, data\$V2, col=data\$V6, xlab = NULL, ylab = NULL, ylim=c(L, 1), xlim=c(1, L), pch = as.numeric(data\$V6))");
	print2file("cmap.R", "\tlegend(\"topright\", bty= \"n\", legend=levels(data\$V6), pch = as.numeric(unique(data\$V6)), col=unique(data\$V6))");
	print2file("cmap.R", "\t} else {");
	print2file("cmap.R", "\t# Plot for a regular RR file");
	print2file("cmap.R", "\tplot(data\$V1, data\$V2, col=data\$V7, xlab = NULL, ylab = NULL, ylim=c(L, 1), xlim=c(1, L), pch = data\$V7)");
	print2file("cmap.R", "\tlegend(\"topright\", bty= \"n\", legend=unique(data\$V8), pch = unique(data\$V7), col = unique(data\$V7))");
	print2file("cmap.R", "}");
	print2file("cmap.R", "text(0, 0, \"$id\", cex=2.0, adj = c(0,1))");
#	print2file("cmap.R", "# display secondary structure information");
#	print2file("cmap.R", "ssfile = \"native.ss\"");
#	print2file("cmap.R", "ss <- file(ssfile, \"rt\")");
#	print2file("cmap.R", "ssline <- readLines(ss, 1)");
#	print2file("cmap.R", "ssvector <- strsplit(ssline, \"\")[[1]]");
#	print2file("cmap.R", "for(i in seq(1,L)){");
#	print2file("cmap.R", "  color = \"cyan\"");
#	print2file("cmap.R", "  if (ssvector[i] == \"H\"){");
#	print2file("cmap.R", "    color = \"orange\"");
#	print2file("cmap.R", "  }");
#	print2file("cmap.R", "  if (ssvector[i] == \"E\"){");
#	print2file("cmap.R", "    color = \"magenta\"");
#	print2file("cmap.R", "  }");
#	print2file("cmap.R", "  arrows(i-1, L+3,  i, L+3, code = 0, lwd = 2, col = color)");
#	print2file("cmap.R", "  arrows( -2, i-1, -2, i, code = 0, lwd = 2, col = color)");
#	print2file("cmap.R", "}");
	
	#print2file("cmap.R", "dev.off()");
	system_cmd("Rscript cmap.R");
}


sub all_pairs_min_dist{
	my $file_pdb = shift;
	my $separation = shift;
	my $d_threshold = shift;
	my $flg_atoms_not_dist = shift;
	confess ":( $file_pdb does not exist!" if not -f $file_pdb;
	my %xyz = ();
	my %pairs_dist = ();
	my %pairs_atoms = ();
	if($param_atom_type eq "ca" or $param_atom_type eq "cb"){
		%xyz = xyz_pdb($file_pdb, $param_atom_type);
		foreach my $r1(sort {$a <=> $b} keys %xyz){
			foreach my $r2(sort {$a <=> $b} keys %xyz){
				next if $r1 >= $r2;
				next if abs($r1 - $r2) < $separation;
				my $d = calc_dist($xyz{$r1}, $xyz{$r2});
				$pairs_dist{$r1." ".$r2} = $d;
				if ($param_atom_type eq "cb"){
					$pairs_atoms{$r1." ".$r2} = "".return_cb_or_ca_atom($r1)." ".return_cb_or_ca_atom($r2);
				}
				else{
					$pairs_atoms{$r1." ".$r2} = "CA CA";
				}
			}
		}
		if ($flg_atoms_not_dist){
			return %pairs_atoms;
		}
		append2log("Returning ".(scalar keys %pairs_dist)." rows for $file_pdb at separation $separation");
		return %pairs_dist;
	}
	else{
		%xyz = xyz_pdb($file_pdb, "ALL");
		foreach my $row1(keys %xyz){
			my @row1 = split /\s+/, $row1;
			my $res1 = $row1[0];
			my $atm1 = $row1[1];
			if ($param_atom_type eq "heavyatoms"){
				next if not ($atm1 eq "N" or $atm1 eq "CA" or $atm1 eq "C" or $atm1 eq "O");  
			}
			foreach my $row2(keys %xyz){
				my @row2 = split /\s+/, $row2;
				my $res2 = $row2[0];
				my $atm2 = $row2[1];
				if ($param_atom_type eq "heavyatoms"){
					next if not ($atm2 eq "N" or $atm2 eq "CA" or $atm2 eq "C" or $atm2 eq "O");  
				}
				next if $res1 >= $res2;
				next if abs($res1 - $res2) < $separation;
				my $d = calc_dist($xyz{$row1}, $xyz{$row2});
				if (not defined $pairs_dist{$res1." ".$res2}){
					$pairs_dist{$res1." ".$res2} = $d;
					$pairs_atoms{$res1." ".$res2} = "$atm1 $atm2" if $flg_atoms_not_dist;
				}
				if ($pairs_dist{$res1." ".$res2} > $d){
					$pairs_dist{$res1." ".$res2} = $d;
					$pairs_atoms{$res1." ".$res2} = "$atm1 $atm2" if $flg_atoms_not_dist;
				}
			}
		}
		if($flg_atoms_not_dist){
			return %pairs_atoms;
		}
		append2log("Returning ".(scalar keys %pairs_dist)." rows for $file_pdb at separation $separation");
		return %pairs_dist;
	}
}

sub append2log{
	my $row = shift;
	print2file($file_log, $row);
}

sub rr_1D_coverage_line{
	my $file_rr = shift;
	my $cont_count = shift;
	error_exit("No RR file!") if not -f $file_rr;
	error_exit("No Contact Count!") if not -f $file_rr;
	my $L = length $rr_sequence;
	my $cov = $rr_sequence;
	$cov =~ s/[A-Z]/-/g;
	my %rr = rr_rows_ordered_in_hash($file_rr, $cont_count, "all", "all");
	foreach (sort {$a <=> $b} keys %rr){
		my @C = split /\s+/, $rr{$_};
		my $r1 = $C[0]; my $r2 = $C[1];
		my $c1 = substr $cov, ($r1 - 1), 1;
		my $c2 = substr $cov, ($r2 - 1), 1;
		if ($c1 eq "-" ){
			$c1 = 1;
		}
		elsif ($c1 eq "*" ){
			$c1 = "*";
		}
		else{
			$c1++;
			$c1 = "*" if ($c1 > 9);
		}
		if ($c2 eq "-" ){
			$c2 = 1;
		}
		elsif ($c2 eq "*" ){
			$c2 = "*";
		}
		else{
			$c2++;
			$c2 = "*" if ($c2 > 9);
		}
		substr $cov, ($r1 - 1), 1, $c1;
		substr $cov, ($r2 - 1), 1, $c2;
	}
	my $cov2 = $cov;
	$cov2 =~ s/-//g;
	return sprintf "$cov [%s pairs %s residues]", (scalar keys %rr), length($cov2);
}

sub print2file{
	my $file = shift;
	my $message = shift;
	my $newline = shift;
	$newline = "\n" if not defined $newline;
	if (-f $file){
		open  FILE, ">>$file" or confess $!;
		print FILE $message."".$newline;
		close FILE;
	}
	else{
		open  FILE, ">$file" or confess $!;
		print FILE $message."".$newline;
		close FILE;
	}
}

sub load_pdb{
	my $dir_chains = shift;
	confess ":( directory $dir_chains does not exist!" if not -d $dir_chains;
	my @pdbList = <$dir_chains/*.pdb>;
	if(not (@pdbList)){
		@pdbList = <$dir_chains/*.ent>;
	}
	confess "ERROR! Directory $dir_chains has no pdb files!\n" unless(@pdbList);
	return @pdbList;
}

sub res_num_res_name{
	my $chain = shift;
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	my %rnum_rname = ();
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		$rnum_rname{parse_pdb_row($_,"rnum")} = parse_pdb_row($_,"rname");
	}
	close CHAIN;
	confess ":(" if not scalar keys %rnum_rname;
	return %rnum_rname;
}

sub pdb2rr{
	my $chain = shift;
	my $rr = shift;
	my $seq_separation = shift;
	my $d_threshold = shift;
	my $param_atom_type = shift;
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	my %contacts_pdb = ();
	my %xyz = xyz_pdb($chain, $param_atom_type);
	confess "\nERROR!! No xyz for any residues in $chain!" if not scalar keys %xyz;
	if (uc($param_atom_type) eq "CA" or uc($param_atom_type) eq "CB"){
		foreach my $r1 (sort keys %xyz){
			foreach my $r2 (sort keys %xyz){
				next if ($r1 >= $r2);
				next if (abs($r2 - $r1) < $seq_separation);
				if ($param_cont_type eq "short"){
					next if (abs($r2 - $r1) > $short_range_max);
				}
				if ($param_cont_type eq "medium"){
					next if (abs($r2 - $r1) > $medium_range_max);
				}
				my $d = calc_dist($xyz{$r1}, $xyz{$r2});
				# Making < $d_threshold and not <= because that is how CASP defines it
				if ($d < $d_threshold){
					$contacts_pdb{"$r1 $r2"} = sprintf "%.3f", $d;
				}
			}
		}
	}
	else{
		foreach my $row1(keys %xyz){
			my @row1 = split /\s+/, $row1;
			my $res1 = $row1[0];
			my $atm1 = $row1[1];
			if ($param_atom_type eq "heavyatoms"){
				next if not ($atm1 eq "N" or $atm1 eq "CA" or $atm1 eq "C" or $atm1 eq "O");  
			}
			foreach my $row2(keys %xyz){
				my @row2 = split /\s+/, $row2;
				my $res2 = $row2[0];
				my $atm2 = $row2[1];
				if ($param_atom_type eq "heavyatoms"){
					next if not ($atm2 eq "N" or $atm2 eq "CA" or $atm2 eq "C" or $atm2 eq "O");  
				}
				next if $res1 >= $res2;
				next if abs($res1 - $res2) < $seq_separation;
				my $d = calc_dist($xyz{$row1}, $xyz{$row2});
				if ($d < $d_threshold){
					$contacts_pdb{"$res1 $res2"} = sprintf "%.3f", $d;
				}
			}
		}
	}
	open RR, ">$rr" or confess $!;
	print RR "".seq_chain_with_gaps($chain)."\n";
	confess "Sorry! There are no such contacts [in $chain] that can be analyzed!</br>Please try assessing ALL instead of LONG-RANGE! That may help!" if not scalar keys %contacts_pdb;
	error_exit("Sorry! There are no such contacts [in $chain] that can be analyzed!</br>Please try assessing ALL instead of LONG-RANGE! That may help!") if not scalar keys %contacts_pdb;
	foreach (sort keys %contacts_pdb){
		print RR "$_ 0 $d_threshold ".$contacts_pdb{$_}."\n";
	}
	close RR;
}
sub seq_chain_with_gaps{
	my $chain = shift;
	my $flag = shift; # flag 1 if the left side dashes of the sequence are not wanted
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	my $start = 1;
	# if flagged, keep start for trimming
	if (defined $flag){
		open CHAIN, $chain or confess $!;
		while(<CHAIN>){
			next if $_ !~ m/^ATOM/;
			if (parse_pdb_row($_,"rname") eq "GLY"){
				next if parse_pdb_row($_,"aname") ne "CA";
			}
			else{
				next if parse_pdb_row($_,"aname") ne "CB";
			}
			confess "ERROR!: ".parse_pdb_row($_,"rname")." residue not defined! \nFile: $chain! \nLine : $_" if (not defined $AA3TO1{parse_pdb_row($_,"rname")});
			$start = parse_pdb_row($_,"rnum");
			last;
		}
		close CHAIN;
	}
	# 1.find end residue number
	my $end;
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		if (parse_pdb_row($_,"rname") eq "GLY"){
			next if parse_pdb_row($_,"aname") ne "CA";
		}
		else{
			next if parse_pdb_row($_,"aname") ne "CB";
		}
		confess "ERROR!: ".parse_pdb_row($_,"rname")." residue not defined! \nFile: $chain! \nLine : $_" if (not defined $AA3TO1{parse_pdb_row($_,"rname")});
		$end = parse_pdb_row($_,"rnum");
	}
	close CHAIN;
	# 2.initialize
	my $seq = "";
	for (my $i = 1; $i <= $end; $i++){
		$seq .= "-";
	}
	# 3.replace with residues
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		if (parse_pdb_row($_,"rname") eq "GLY"){
			next if parse_pdb_row($_,"aname") ne "CA";
		}
		else{
			next if parse_pdb_row($_,"aname") ne "CB";
		}
		confess "ERROR!: ".parse_pdb_row($_,"rname")." residue not defined! \nFile: $chain! \nLine : $_" if (not defined $AA3TO1{parse_pdb_row($_,"rname")});
		my $rnum = parse_pdb_row($_,"rnum");
		$rnum =~ s/[A-G]//g;
		substr $seq, ($rnum - 1), 1, $AA3TO1{parse_pdb_row($_,"rname")}; 
	}
	close CHAIN;
	confess "$chain has less than 1 residue!" if (length($seq) < 1);
	return (substr $seq, $start - 1);
}
sub seq_chain{
	my $chain = shift;
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	my $seq = "";
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		if (parse_pdb_row($_,"rname") eq "GLY"){
			next if parse_pdb_row($_,"aname") ne "CA";
		}
		else{
			next if parse_pdb_row($_,"aname") ne "CB";
		}
		confess "ERROR!: ".parse_pdb_row($_,"rname")." residue not defined! \nFile: $chain! \nLine : $_" if (not defined $AA3TO1{parse_pdb_row($_,"rname")});
		my $res = $AA3TO1{parse_pdb_row($_,"rname")};
		$seq .= $res;
	}
	close CHAIN;
	confess "$chain has less than 1 residue!" if (length($seq) < 1);
	return $seq;
}
sub parse_pdb_row{
	my $row = shift;
	my $param = shift;
	my $result;
	$result = substr($row,6,5) if ($param eq "anum");
	$result = substr($row,12,4) if ($param eq "aname");
	$result = substr($row,16,1) if ($param eq "altloc");
	$result = substr($row,17,3) if ($param eq "rname");
	$result = substr($row,22,5) if ($param eq "rnum");
	$result = substr($row,26,1) if ($param eq "insertion");
	$result = substr($row,21,1) if ($param eq "chain");
	$result = substr($row,30,8) if ($param eq "x");
	$result = substr($row,38,8) if ($param eq "y");
	$result = substr($row,46,8) if ($param eq "z");
	confess "Invalid row[$row] or parameter[$param]" if (not defined $result);
	$result =~ s/\s+//g;
	return $result;
}

sub reindex_chain{
	my $file_pdb = shift;
	my $index = shift;
	my $out_pdb = shift;
	confess "ERROR! file $file_pdb does not exist!" if not -f $file_pdb;
	confess "ERROR! index $index is invalied!" if not defined $index;
	open PDBFILE, $file_pdb or confess $!;
	my @pdb_lines = <PDBFILE>;
	close PDBFILE;
	# (c) Reindex Chain. Assumptions: non-standard residues removed, alternative locations removed, one model, one chain.
	my $res_counter = $index - 1;
	my $atom_counter = 0;
	my $prev_res_num = "XX";
	open OUTPDB, ">$out_pdb" or confess $!;
	foreach (@pdb_lines) {
		next if $_ !~ m/^ATOM/;
		next if not ((parse_pdb_row($_,"altloc") eq "") or (parse_pdb_row($_,"altloc") eq "A"));
		next if not defined $AA3TO1{parse_pdb_row($_,"rname")};
		my $this_rnum = parse_pdb_row($_,"rnum");
		if ($prev_res_num ne $this_rnum) {
			$prev_res_num = $this_rnum;
			$res_counter++;
		}
		$atom_counter++;
		my $rnum_string = sprintf("%4s", $res_counter);
		my $anum_string = sprintf("%5s", $atom_counter);
		my $row = substr($_,0,6).$anum_string.substr($_,11,5)." ".substr($_,17,3)." "." ".$rnum_string." ".substr($_,27);
		print OUTPDB $row;
	}
	print OUTPDB "END\n";
	close OUTPDB;
}

sub calc_dist{
	my $x1y1z1 = shift;
	my $x2y2z2 = shift;
	my @row1 = split(/\s+/, $x1y1z1);
	my $x1 = $row1[0]; my $y1 = $row1[1]; my $z1 = $row1[2];
	my @row2 = split(/\s+/, $x2y2z2);
	my $x2 = $row2[0]; my $y2 = $row2[1]; my $z2 = $row2[2];
	my $d = sprintf "%.3f", sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
	return $d;
}
sub rrfile_to_r1r2hash{
	my $file_rr = shift;
	my $seq_sep = shift;
	my $count = shift;
	$seq_sep = 1 if not defined $seq_sep;
	$count = 1000000 if not defined $count;
	confess ":(" if not -f $file_rr;
	my %contacts = ();
	open RR, $file_rr or confess $!;
	while(<RR>){
		next unless ($_ =~ /[0-9]/);
		$_ =~ s/^\s+//;
		next unless ($_ =~ /^[0-9]/);
		my @C = split(/\s+/, $_);
		confess "ERROR! Expecting a pair in row [".$_."]!\n" if (not defined $C[0] || not defined $C[1]);
		next if (abs($C[1] - $C[0]) < $seq_sep);
		if(defined $C[4]){
			$contacts{$C[0]." ".$C[1]} = $C[4];
		}
		elsif(defined $C[2] && $C[2] != 0){
			$contacts{$C[0]." ".$C[1]} = $C[2];
		}
		else{
			confess "ERROR! Confidence column not defined in row [".$_."] in file $file_rr!\n";
		}
		last if (scalar keys %contacts) == $count;
	}
	close RR;
	return %contacts;
}

sub rr_has_seq{
	my $file_rr = shift;
	confess "ERROR! Input file $file_rr does not exist!" if not -f $file_rr;
	my $seq;
	open RR, $file_rr or confess "ERROR! Could not open $file_rr! $!";
	while(<RR>){
		chomp $_;
		$_ =~ tr/\r//d; # chomp does not remove \r
		$_ =~ s/^\s+//;
		next if ($_ =~ /^PFRMAT/);
		next if ($_ =~ /^TARGET/);
		next if ($_ =~ /^AUTHOR/);
		next if ($_ =~ /^SCORE/); 
		next if ($_ =~ /^REMARK/);
		next if ($_ =~ /^METHOD/);
		next if ($_ =~ /^MODEL/); 
		next if ($_ =~ /^PARENT/);
		last if ($_ =~ /^TER/);   
		last if ($_ =~ /^END/);
		# Now, I can directly merge to RR files with sequences on top
		last if ($_ =~ /^[0-9]/);
		$seq .= $_;
	}
	close RR;
	return 0 if not defined $seq;
	return 1;
}

sub seq_rr{
	my $file_rr = shift;
	confess "ERROR! Input file $file_rr does not exist!" if not -f $file_rr;
	my $seq;
	open RR, $file_rr or confess "ERROR! Could not open $file_rr! $!";
	while(<RR>){
		chomp $_;
		$_ =~ s/\r//g; # chomp does not remove \r
		$_ =~ s/^\s+//;
		$_ =~ s/\s+//g;
		next if ($_ =~ /^>/);
		next if ($_ =~ /^PFRMAT/);
		next if ($_ =~ /^TARGET/);
		next if ($_ =~ /^AUTHOR/);
		next if ($_ =~ /^SCORE/); 
		next if ($_ =~ /^REMARK/);
		next if ($_ =~ /^METHOD/);
		next if ($_ =~ /^MODEL/); 
		next if ($_ =~ /^PARENT/);
		last if ($_ =~ /^TER/);   
		last if ($_ =~ /^END/);
		# Now, I can directly merge to RR files with sequences on top
		last if ($_ =~ /^[0-9]/);
		$seq .= $_;
	}
	close RR;
	error_exit("Input RR file does not have sequence row!</br>RR-file : <b>$file_rr</b></br>Please make sure that all input RR files have sequence headers!") if not defined $seq;
	return $seq;
}

sub ss2residues_hash{
	my $file_ss = shift;
	confess ":(" if not -f $file_ss;
	my $seq = seq_ss($file_ss);
	my %res = ();
	foreach (my $i = 0; $i <= length($seq); $i++){
		$res{$i+1} = substr $seq, $i, 1;
	}
	return %res;
}
sub seq_fasta{
	my $file_fasta = shift;
	confess ":(" if not -f $file_fasta;
	my $seq = "";
	open FASTA, $file_fasta or confess $!;
	while (<FASTA>){
		next if (substr($_,0,1) eq ">"); 
		chomp $_;
		$_ =~ tr/\r//d; # chomp does not remove \r
		$seq .= $_;
	}
	close FASTA;
	return $seq;
}

sub fasta2residues_hash{
	my $file_fasta = shift;
	confess ":(" if not -f $file_fasta;
	my $seq = seq_fasta($file_fasta);
	my %res = ();
	foreach (my $i = 0; $i < length($seq); $i++){
		$res{$i+1} = substr $seq, $i, 1;
	}
	return %res;
}

sub return_cb_or_ca_atom{
	my $rnum = shift;
	my $residue = substr $rr_sequence, $rnum - 1, 1;
	confess "rnum not defined!" if not defined $rnum;
	confess "Could not find residue name for $rnum!" if not $residue;
	return "CA" if $residue eq "G";
	return "CB";
}

sub xyz_pdb{
	my $chain = shift;
	my $atom_selection = shift; # ca or cb or all/any
	$atom_selection = "all" if $atom_selection eq "any";
	confess "\nERROR! file $chain does not exist!" if not -f $chain;
	confess "\nERROR! Selection must be ca or cb or all or heavyatoms" if not (uc($atom_selection) eq "CA" or uc($atom_selection) eq "ALL" or uc($atom_selection) eq "CB" or uc($atom_selection) eq "HEAVYATOMS");
	my %xyz_pdb = ();
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		$xyz_pdb{"".parse_pdb_row($_,"rnum")." ".parse_pdb_row($_,"aname")} = "".parse_pdb_row($_,"x")." ".parse_pdb_row($_,"y")." ".parse_pdb_row($_,"z");
	}
	close CHAIN;
	confess "\nERROR!: xyz_pdb is empty\n" if (not scalar keys %xyz_pdb);
	if (uc($atom_selection) eq "ALL"){
		return %xyz_pdb;
	}
	elsif (uc($atom_selection) eq "HEAVYATOMS"){
		foreach (sort keys %xyz_pdb){
			my @C = split /\s+/, $_;
			if (not($C[1] eq "N" or $C[1] eq "CA" or $C[1] eq "C" or $C[1] eq "O")){
				delete $xyz_pdb{$_};
			}
		}
		return %xyz_pdb;
	}
	my %native_res_list = res_num_res_name($chain);
	my %selected_xyz = ();
	foreach (sort keys %xyz_pdb){
		my @C = split /\s+/, $_;
		my $atom_of_interest = uc($atom_selection);
		$atom_of_interest = "CA" if $native_res_list{$C[0]} eq "GLY";
		# Some pdb files have errors. Non-gly residues do not have CB atoms.
		# For instance, http://sysbio.rnet.missouri.edu/coneva/preloaded_data/fragfold/native/1ej8A_reindexed.pdb
		# Need to throw errors in such cases.
		error_exit("The pdb file $param_pdb_file does not have CB atom for residue ".$C[0]."! Try assessing CA-contacts instead!") if not defined $xyz_pdb{$C[0]." ".$atom_of_interest};
		next if $C[1] ne $atom_of_interest;
		$selected_xyz{$C[0]} = $xyz_pdb{$_};
	}
	confess "\nERROR! Empty xyz coordinates in the pdb file!" if not scalar keys %selected_xyz;
	return %selected_xyz;
}

sub system_cmd{
	my $command = shift;
	error_exit("EXECUTE [$command]") if (length($command) < 5  and $command =~ m/^rm/);
	append2log("[Executing $command]");
	system($command);
	if($? != 0){
		my $exit_code  = $? >> 8;
		error_exit("Failed executing [$command]!<br>ERROR: $!");
	}
}

sub set_progress{
	my $percent = shift;
	print "\n<script type=\"text/javascript\">updateProgress(\"$percent\");</script>";
}
