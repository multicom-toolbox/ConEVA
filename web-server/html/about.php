<?php include("header.php");?>
<table id="main" frame="box" width="1000" cellspacing="5" cellpadding="5" border="0" rules="none" align="center" bgcolor="#F0F0F0 ">
<tr><td align="center">
		<table style="display:" id="usage" class="results">
		<th align=left><font style="color:black;font-size:150%">Usage</font><th>
		<tr><td style="text-align:left;">
			<ul><li>
			<b>CASP RR format</b> description is available at <a href="http://predictioncenter.org/casprol/index.cgi?page=format#RR">CASP's website<a>.
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
			When a zip file of RR files is submitted as input, CONEVA can handle maximum 10 RR files (in the current settings). If the submitted zip file has more than 10 files, top 10 (sorted alphabetically) files are only processed.
			</li><li>
			CONEVA computes results within 2 minutes for most of the inputs. Longer proteins and computationally intensive contact definitions (like all-atom definitions) take longer time.
			</li>
			</ul>
		</td></tr>
		<th align=left><font style="color:black;font-size:150%">Screenshots</font><th>
		<tr><td style="text-align:left;">
			<ul>
			<li><a href="http://cactus.rnet.missouri.edu/coneva/screenshots/">Images of input and output sample</a>
			</li>
			</ul>
		</td></tr>

	</table>
<!--
	<b>Authors:</b></br>
	<a href="https://sites.google.com/site/visitbadri/">Badri Adhikari</a></br>
	Debswapna Bhattacharya</br>
	Renzhi Cao</br>
	Jie Hou</br>
	<a href="http://calla.rnet.missouri.edu/cheng/">Dr. Jianlin Cheng (PI)</a></br>
	</br>
	<b>Paper:</b></br>
	<a href="" style="text-decoration: none">CONASSESS: an interactive web server for a comprehensive assessment of predicted protein contacts</a>,</br> 
	(In preparation)</br>
	B. Adhikari, D. Bhattacharya, J.Hou, J. Cheng. </br>
	</br>
	<b>Book Chapter:</b></br>
	<a href="" style="text-decoration: none">Assessing predicted contacts for building protein three-dimensional models </a>,</br> 
	(In press)</br>
	B. Adhikari, D. Bhattacharya, R. Cao, J. Cheng. </br>
-->	</br>
	Bioinformatics and Systems Biology Laboratory (<a href="http://calla.rnet.missouri.edu/cheng/cheng_research.html">BDM</a>)</a></br>
	Department of Computer Science</a></br>   
	University of Missouri</a></br> Columbia, MO</br>
	</p>
</td></tr>
<?php include("footer.php");?>
