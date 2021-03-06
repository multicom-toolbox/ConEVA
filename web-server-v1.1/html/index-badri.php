<?php include("header.php");?>
	<form method="post" onsubmit="return validateForm()" action="http://cactus.rnet.missouri.edu/cgi-bin/conassess/main_v2.0.cgi" enctype="multipart/form-data">
	<table id="main" frame="box" width="1000" cellspacing="5" cellpadding="5" border="0" rules="none" align="center" bgcolor="#F0F0F0">
		<tr id="tr_rr_raw">
			<td valign="top" align="right">Contacts</br>(or Zip of RRs)</td>
			<td><textarea id="rr_raw" name="rr_raw" rows="6" cols="100" placeholder="(1)Paste RR file contents or 
(2)Paste a http link to a RR file"></textarea></td>
		</tr>
		<tr id="tr_rr_upload">
			<td></td>
			<td>Upload a RR file or a Zip file <input type="file" id="rrfile" name="rrfile" size="100"><hr></td>
		</tr>
		<tr>
			<td valign="top" align="right">PDB</td>
			<td><textarea id="pdb_raw" name="pdb_raw" rows="4" cols="100"  placeholder="(1)Insert PDB id using 'id:pdb-id', for instance ‘id:3e7u’, or
(2)Paste a http link, for instance, 'http://www.rcsb.org/pdb/files/3E7U.pdb' or
(3)Paste contents of PDB file"></textarea></td>
		</tr>
		<tr>
			<td></td>
			<td><font>Upload a PDB structure file <input type="file" id="pdbfile" name="pdbfile" size="100"></font></td>
		</tr>
		<tr>
			<td align="right"></td>
			<td><input id="native_rr_flag" onchange="check_and_show_dthres_row()" type="checkbox" name="native_rr_flag"><font >"Analyze PDB's contacts" calculating contacts from the PDB (ignores RR input)</font><hr></td>
		</tr>
		<tr id="distthresholdrow" style="display:none">
			<td valign="top" align="right">Distance Threshold</td>
			<td style="font-family:monospace">
			<select id="distthreshold" name="distthreshold">
				<option value="6">6</option>
				<option value="7">7</option>
				<option value="8" selected>8</option>
				<option value="9">9</option>
				<option value="10">10</option>
				<option value="11">11</option>
				<option value="12">12</option>
			</select>&nbsp;&#197;
			</td>
		</tr>	
	<!--	<tr>
			<td valign="top" align="right">Secondary Structure</td>
			<td><textarea id="protein_sec" name="protein_sec" rows="3" cols="100" placeholder='(optional) Paste 3 state SS sequence like CCCHHHHCCCCCEEEECCCCEEEECCC'></textarea></td>
		</tr>
	-->	<tr>
			<td valign="top" align="right">Contact Atoms</td>
			<td><font face="monospace">
			<input type="radio" name="atomtype" id="atomtypecb"value="cb" checked>C&beta;
			<input type="radio" name="atomtype" id="atomtypeca" value="ca">C&alpha;
			<input type="radio" name="atomtype" id="atomtypeheavy" value="heavyatoms">N/C/C&alpha;/O (any backbone atoms)
			<input type="radio" name="atomtype" id="atomtypeany" value="any">any (very slow)
			</font>
			</td>
		</tr>
		<tr>
			<td valign="top" align="right">Contact Type</td>
			<td><font face="monospace">
			<input type="radio" name="contacttype" id="contacttypeall"    value="all" checked>All
			<input type="radio" name="contacttype" id="contacttypelong"   value="long">Long-Range
			<input type="radio" name="contacttype" id="contacttypemedium" value="medium">Medium-Range
			<input type="radio" name="contacttype" id="contacttypeshort"  value="short">Short-Range
			</td>
		</tr>
		<tr>
			<td valign="top" align="right">Sequence Separation</td>
				<td> <font face="monospace">
				Short-range &nbsp;min <input id="srmin" name="srmin" size="5" value=6 style="font-family:monospace" required/> max <input id="srmax" name="srmax" size="5" value=11 style="font-family:monospace" required/> </br>
				Medium-range min <input id="mrmin" name="mrmin" size="5" value=12 style="font-family:monospace" required/> max <input id="mrmax" name="mrmax" size="5" value=23 style="font-family:monospace" required/> </br>
				Long-range &nbsp;&nbsp;min <input id="lrmin" name="lrmin" size="5" value=24 style="font-family:monospace" required/> max <input id="lrmax" name="lrmax" size="5" value=10000 style="font-family:monospace" required/> </br>
				</font>
			</td>
		</tr>
	<!--
		<tr>
		<td valign="top" align="right">Max. No. of Top Ranked Contacts to load from each RR file</td>
			<td> <font face="monospace">
				<input type="radio" name="maxcon" id="maxconL"   value="L">L
				<input type="radio" name="maxcon" id="maxcon2L"  value="2L">2L
				<input type="radio" name="maxcon" id="maxconALL" value="ALL" checked>ALL (slower)
				</br>
				</font>
			</td>
		</tr>
	-->
		<tr>
		<td valign="top" align="right">Neighbor relaxation for </br>Jaccard similarity calculations</td>
			<td> <font face="monospace">
				<input type="radio" name="neighborsize" id="neighborsize0" value="0" checked>0
				<input type="radio" name="neighborsize" id="neighborsize1" value="1">1
				<input type="radio" name="neighborsize" id="neighborsize2" value="2">2
				<input type="radio" name="neighborsize" id="neighborsize3" value="3">3
				</br>
				</font>
			</td>
		</tr>
		<tr><td colspan=2 align=center>
			<input id="submitbutton" name="action" type="submit" style="font-size:18px" value="Assess"/>
		</td></tr>
		<tr>
			<td></td><td style="text-align:right; color:blue;"><u><a onclick="toggle_it('examples');toggle_it('usage')">Pre-curated Data sets</a></u> &nbsp; &nbsp;<u><a onclick="clearForm()">Reset Fields</a></u></td>
		</tr>
	</table>
	<font>
	<table style="display:" id="usage" class="results">
		<th align=left><font style="color:black;font-size:150%">Usage</font><th>
		<tr><td style="text-align:left;">
			<ul><li>
			<b>CASP RR format</b> description is available at <a href="http://predictioncenter.org/casprol/index.cgi?page=format#RR">CASP's website<a>.
			</li><li>
			<b>Native PDB</b> may be supplied in 4 ways: (1) a four or five letter PDB id using 'id:pdbid' (b) an http link of the pdb file (c) raw pdb content pasted in the text field (d) Upload a pdb file in local computer.
			</li><li>
			Checking the '<b>Analyze PDB's contacts</b>' option takes the input PDB file, calculates contacts (based on the definition selected), and evaluates the 'true' contacts.
			</li><li>
			<b>Contact atoms</b> defines the atoms used to define contacts in RR file. For instance, CASP uses Cβ atom for defining contacts.
			</li><li>
			<b>Contact type</b> defines the type of contacts to be selected for analysis; other contacts are ignored. For instance, CASP focuses on Long-Range contact analysis.
			</li><li>
			<b>Sequence Separation</b> sets minimum and maximum sequence separation distances for defining short-, medium-, and long-range contacts. For instance, for analyzing CASP Long-Range contacts, set the long-range min to 24 and max to a very high value, and select Long-Range as contact type.
			</li><li>
			<b>Neighbor relaxation for Jaccard similarity coefficient</b> sets neighborhood size for relaxing the contact definition. It is used only when more than one set of contacts are supplied as input. If it is set to 0 (default) contact defintions are not changed. If it is set to 1, contacts with residue separation ±1 are considered same while computing set unions and intersections. Higher values will find deeper similarity between two sets of contacts.
			</ul>
		</td></tr>
		<th align=left><font style="color:black;font-size:150%">Notes</font><th>
		<tr><td style="text-align:left;">
			<ul>
			<li>
			Recommended screen resolution is 1280 X 720 or higher.
			</li><li>
			When a zip file of RR files is submitted as input, CONASSESS can handle maximum 10 RR files (in the current settings). If the submitted zip file has more than 10 files, top 10 (sorted alphabetically) files are only processed.
			</li><li>
			CONASSESS computes results within 2 minutes for most of the inputs. Longer proteins and computationally intensive contact definitions (like all-atom definitions) take longer time.
			</li>
			</ul>
		</td></tr>
	</table>
	<font style="color:blue;">
	 <table style="display:none" id="examples" class="results">
		<tr><td></td></tr>
		<tr><td></td></tr>
		<th align=left><font style="color:black;font-size:120%">FRAGFOLD data set of 150 proteins</font><th>
		<tr><td style="text-align:left;">
			<font style="color: black;"><b>Set 1: Evaluate and Compare contacts in MetaPSICOV <a href="http://bioinfadmin.cs.ucl.ac.uk/downloads/MetaPSICOV/" target="_blank">supplementary data<a>. DNcon predictions were made locally.</br></b></font>
			<?php
				error_reporting(E_ALL);
				ini_set('display_errors', 1);
				$handle = @fopen("preloaded/lists/fragfold.txt", "r");
				if ($handle) {
					while (($buffer = fgets($handle, 4096)) !== false) {
						$id = rtrim($buffer);
						$rrurl = "http://cactus.rnet.missouri.edu/conassess/preloaded/metapsicov/merged_n_zipped/$id.zip";
						$nativeurl = "http://cactus.rnet.missouri.edu/conassess/preloaded/metapsicov/pdb_frg_supp/$id.pdb";
						echo "<a style=\"text-decoration:none;\" onclick=\"evaluate_rr('metapsicov', '$id', '$rrurl', '$nativeurl')\">$id</a> \n";
					}
					if (!feof($handle)) {
						echo "Error: unexpected fgets() fail\n";
					}
					fclose($handle);
				}
			?> 
			</td></tr>
		<tr><td style="text-align:left;">
			<font style="color: black;"><b>Set 2: Analyze the native PDB's contacts.</br></b></font>
			<?php
				error_reporting(E_ALL);
				ini_set('display_errors', 1);
				$handle = @fopen("preloaded/lists/fragfold.txt", "r");
				if ($handle) {
					while (($buffer = fgets($handle, 4096)) !== false) {
						$id = rtrim($buffer);
						echo "<a style=\"text-decoration:none;\" onclick=\"native_analysis('metapsicov', '$id', 'id:$id')\">$id</a> \n";
					}
					if (!feof($handle)) {
						echo "Error: unexpected fgets() fail\n";
					}
					fclose($handle);
				}
			?>
			</td> 
		</tr>
		<tr><td style="text-align:left;">
			<font style="color: black;"><b>Set 3: Evaluate contacts predicted by PSICOV only.</br></b></font>
			<?php
				error_reporting(E_ALL);
				ini_set('display_errors', 1);
				$handle = @fopen("preloaded/lists/fragfold.txt", "r");
				if ($handle) {
					while (($buffer = fgets($handle, 4096)) !== false) {
						$id = rtrim($buffer);
						$rrurl = "http://cactus.rnet.missouri.edu/conassess/preloaded/metapsicov/casp_format/psicov/$id.contacts";
						$nativeurl = "http://cactus.rnet.missouri.edu/conassess/preloaded/metapsicov/pdb_frg_supp/$id.pdb";
						echo "<a style=\"text-decoration:none;\" onclick=\"evaluate_rr('psicov', '$id', '$rrurl', '$nativeurl')\">$id</a> \n";
					}
					if (!feof($handle)) {
						echo "Error: unexpected fgets() fail\n";
					}
					fclose($handle);
				}
			?>
			</td></tr>
		<tr><td></td></tr>
		<tr><td></td></tr>
		<tr><td></td></tr>
		<th align=left><font style="color:black;font-size:120%">Some proteins in EVFOLD data set of 15 proteins</font><th>
		<tr><td style="text-align:left;">
			<font style="color: black;"><b>Set 4: Evaluate contacts predicted by EVFOLD method (<a href="http://evfold.org/foldingproteins/evfold/download/StructureFromSequences_Appendix_A1.tar.gz" target="_blank">supp</a>).</br></b></font>
			<?php
				error_reporting(E_ALL);
				ini_set('display_errors', 1);
				$handle = @fopen("preloaded/lists/evfold.txt", "r");
				if ($handle) {
					while (($buffer = fgets($handle, 4096)) !== false) {
						$pair = rtrim($buffer);
						$newpair = explode(" ", $pair, 2);
						$rrurl = "http://cactus.rnet.missouri.edu/conassess/preloaded/evfold/input/$newpair[0].rr";
						$nativeurl = "http://cactus.rnet.missouri.edu/conassess/preloaded/evfold/native/$newpair[1]_filtered.pdb";
						$prompt = $newpair[0]."-".$newpair[1];
						echo "<a style=\"text-decoration:none;\" onclick=\"evaluate_rr('evfold', '$prompt', '$rrurl', '$nativeurl')\">$newpair[0]</a> &nbsp;&nbsp;\n";
					}
					if (!feof($handle)) {
						echo "Error: unexpected fgets() fail\n";
					}
					fclose($handle);
				}
			?>
		</td></tr>
		<tr><td></td></tr>
		<tr><td></td></tr>
		<tr><td></td></tr>
		<th align=left><font style="color:black;font-size:120%">CASP data sets</font><th>
		<tr><td style="text-align:left;">
			<font style="color: black;"><b>Set 5: CASP11 RR - Evaluate and compare contacts predicted by various contact prediction groups.</br></b></font>
			<?php
				error_reporting(E_ALL);
				ini_set('display_errors', 1);
				$handle = @fopen("preloaded/lists/casp11_targets.txt", "r");
				if ($handle) {
					while (($buffer = fgets($handle, 4096)) !== false) {
						$id = rtrim($buffer);
						$rrurl = "http://cactus.rnet.missouri.edu/conassess/preloaded/casp/11/RR/predictions/$id/$id.zip";
						$nativeurl = "http://cactus.rnet.missouri.edu/conassess/preloaded/casp/11/RR/native/$id.pdb";
						echo "<a style=\"text-decoration:none;\" onclick=\"evaluate_rr('casprr', '$id', '$rrurl', '$nativeurl')\">$id</a> \n";
					}
					if (!feof($handle)) {
						echo "Error: unexpected fgets() fail\n";
					}
					fclose($handle);
				}
			?>
			</td></tr>
		<tr><td style="text-align:left;">
			<font style="color: black;"><b>Set 6: Evaluate contacts predicted by the group 038 (can be changed to any group of interest) in CASP11 RR category.</br></b></font>
			<?php
				error_reporting(E_ALL);
				ini_set('display_errors', 1);
				$handle = @fopen("preloaded/lists/casp11_targets.txt", "r");
				if ($handle) {
					while (($buffer = fgets($handle, 4096)) !== false) {
						$id = rtrim($buffer);
						$rrurl = "http://cactus.rnet.missouri.edu/conassess/preloaded/casp/11/RR/predictions/$id/${id}RR038_1";
						$nativeurl = "http://cactus.rnet.missouri.edu/conassess/preloaded/casp/11/RR/native/$id.pdb";
						echo "<a style=\"text-decoration:none;\" onclick=\"evaluate_rr('casprr', '$id by group 038', '$rrurl', '$nativeurl')\">$id</a> \n";
					}
					if (!feof($handle)) {
						echo "Error: unexpected fgets() fail\n";
					}
					fclose($handle);
				}
			?>
		</td></tr>
		<tr><td style="text-align:left;">
			<font style="color: black;"><b>Set 7: CASP11 RR - Evaluate and compare contacts predicted by selected contact prediction groups.</br></b></font>
			<?php
				error_reporting(E_ALL);
				ini_set('display_errors', 1);
				$handle = @fopen("preloaded/lists/casp11_targets.txt", "r");
				if ($handle) {
					while (($buffer = fgets($handle, 4096)) !== false) {
						$id = rtrim($buffer);
						$rrurl = "http://cactus.rnet.missouri.edu/conassess/preloaded/casp/11/RR_selected/zipped/$id.zip";
						$nativeurl = "http://cactus.rnet.missouri.edu/conassess/preloaded/casp/11/RR/native/$id.pdb";
						echo "<a style=\"text-decoration:none;\" onclick=\"evaluate_rr('casprr', '$id', '$rrurl', '$nativeurl')\">$id</a> \n";
					}
					if (!feof($handle)) {
						echo "Error: unexpected fgets() fail\n";
					}
					fclose($handle);
				}
			?>
			</td></tr>
		<tr><td style="text-align:left;">
			<font style="color: black;"><b>Set 8: Analyze the native PDB domains of CASP10.</br></b></font>
			<?php
				error_reporting(E_ALL);
				ini_set('display_errors', 1);
				$handle = @fopen("preloaded/lists/casp10_domains.txt", "r");
				if ($handle) {
					while (($buffer = fgets($handle, 4096)) !== false) {
						$id = rtrim($buffer);
						$nativeurl = "http://cactus.rnet.missouri.edu/conassess/preloaded/casp/10/native_domains/$id.pdb";
						echo "<a style=\"text-decoration:none;\" onclick=\"native_analysis('casp', '$id', '$nativeurl')\">$id</a> \n";
					}
					if (!feof($handle)) {
						echo "Error: unexpected fgets() fail\n";
					}
					fclose($handle);
				}
			?>
			</td></tr>
		<tr><td style="text-align:left;">
			<font style="color: black;"><b>Set 9: Analyze the native PDB domains of CASP11.</br></b></font>
			<?php
				error_reporting(E_ALL);
				ini_set('display_errors', 1);
				$handle = @fopen("preloaded/lists/casp11_domains.txt", "r");
				if ($handle) {
					while (($buffer = fgets($handle, 4096)) !== false) {
						$id = rtrim($buffer);
						$nativeurl = "http://cactus.rnet.missouri.edu/conassess/preloaded/casp/11/native_domains/$id.pdb";
						echo "<a style=\"text-decoration:none;\" onclick=\"native_analysis('casp', '$id', '$nativeurl')\">$id</a> \n";
					}
					if (!feof($handle)) {
						echo "Error: unexpected fgets() fail\n";
					}
					fclose($handle);
				}
			?>
			</td></tr>
			<tr><td style="text-align:left;">
			<font style="color: black;"><b>Set 10: Compare DNcon and metaPSICOV on 15 CASP11 FM targets (metaPSICOV was run locally).</br></b></font>
			<?php
				error_reporting(E_ALL);
				ini_set('display_errors', 1);
				$handle = @fopen("preloaded/lists/casp11_fm_15.txt", "r");
				if ($handle) {
					while (($buffer = fgets($handle, 4096)) !== false) {
						$id = rtrim($buffer);
						$rrurl = "http://cactus.rnet.missouri.edu/conassess/preloaded/dncon_vs_metaPSICOV/$id.zip";
						$nativeurl = "http://cactus.rnet.missouri.edu/conassess/preloaded/casp/11/RR/native/$id.pdb";
						echo "<a style=\"text-decoration:none;\" onclick=\"evaluate_rr('casprr', '$id', '$rrurl', '$nativeurl')\">$id</a> \n";
					}
					if (!feof($handle)) {
						echo "Error: unexpected fgets() fail\n";
					}
					fclose($handle);
				}
			?>
			</td></tr>
			<tr><td style="text-align:left;">
			<font style="color: black;"><b>Set 11: GDFuzz3D - 45 CASP10 TBM targets.</br></b></font>
			<?php
				error_reporting(E_ALL);
				ini_set('display_errors', 1);
				$handle = @fopen("preloaded/lists/GDFuzz3D_casp10_tbm45.txt", "r");
				if ($handle) {
					while (($buffer = fgets($handle, 4096)) !== false) {
						$id = rtrim($buffer);
						$rrurl = "http://cactus.rnet.missouri.edu/conassess/preloaded/GDFuzz3D/multicom_fuzzy_cm/$id.txt";
						$nativeurl = "http://cactus.rnet.missouri.edu/conassess/preloaded/casp/10/native_targets/$id.pdb";
						echo "<a style=\"text-decoration:none;\" onclick=\"evaluate_rr('gdfuzz3d', '$id', '$rrurl', '$nativeurl')\">$id</a> \n";
					}
					if (!feof($handle)) {
						echo "Error: unexpected fgets() fail\n";
					}
					fclose($handle);
				}
			?>
			</td></tr>
			<tr><td style="text-align:left;">
			<font style="color: black;"><b>Set 12: GDFuzz3D - 10 CASP10 TBM targets (with good contact maps).</br></b></font>
			<?php
				error_reporting(E_ALL);
				ini_set('display_errors', 1);
				$handle = @fopen("preloaded/lists/GDFuzz3D_casp10_tbm10.txt", "r");
				if ($handle) {
					while (($buffer = fgets($handle, 4096)) !== false) {
						$id = rtrim($buffer);
						$rrurl = "http://cactus.rnet.missouri.edu/conassess/preloaded/GDFuzz3D/multicom_fuzzy_cm/$id.txt";
						$nativeurl = "http://cactus.rnet.missouri.edu/conassess/preloaded/casp/10/native_targets/$id.pdb";
						echo "<a style=\"text-decoration:none;\" onclick=\"evaluate_rr('gdfuzz3d', '$id', '$rrurl', '$nativeurl')\">$id</a> \n";
					}
					if (!feof($handle)) {
						echo "Error: unexpected fgets() fail\n";
					}
					fclose($handle);
				}
			?>
			</td></tr>
		</table>
	</form>
<?php include("footer.php");?>