function clearForm(){
	document.getElementById("rr_raw").value            = "";
	document.getElementById("pdb_raw").value           = "";
	document.getElementById("rrfile").value            = "";
	document.getElementById("pdbfile").value           = "";
	document.getElementById("native_rr_flag").checked  = false;
	document.getElementById("srmin").value             = 6;
	document.getElementById("mrmin").value             = 12;
	document.getElementById("lrmin").value             = 24;
	document.getElementById("srmax").value             = 11;
	document.getElementById("mrmax").value             = 23;
	document.getElementById("lrmax").value             = 10000;
	document.getElementById("contacttypeall").checked  = true;
	document.getElementById("atomtypecb").checked      = true;
	document.getElementById("neighborsize0").checked   = true;
	document.getElementById('distthresholdrow').style.display = 'none';
//	check_and_show_dthres_row();
}

window.onload = function () {
	check_and_show_dthres_row();
}

function evaluate_rr(data, prompt, rrurl, nativeurl){
	var r = confirm("Load data for ".concat(prompt, "?"));
	if (r == false) {
		return;
	}
	scrollToTop(500);
	clearForm();
	document.getElementById("rr_raw").value = rrurl;
	document.getElementById("pdb_raw").value = nativeurl;
	if(data == "metapsicov"){
		document.getElementById("srmin").value = 5;
		document.getElementById("contacttypelong").checked = true;
	}
	if(data == "casprr"){
		document.getElementById("srmin").value = 6;
		document.getElementById("contacttypelong").checked = true;
	}
	if(data == "gdfuzz3d"){
		document.getElementById("srmin").value = 6;
		document.getElementById("contacttypeany").checked = true;
	}
	if(data == "psicov"){
		document.getElementById("srmin").value = 5;
		document.getElementById("contacttypeall").checked = true;
	}
	if(data == "evfold"){
		document.getElementById("srmin").value = 5;
		document.getElementById("contacttypeall").checked = true;
		document.getElementById("atomtypeany").checked     = true;
	}
	document.getElementById('submitbutton').focus();
}

function assess_rr(data, prompt, rrurl){
	var r = confirm("Load data for ".concat(prompt, "?"));
	if (r == false) {
		return;
	}
	scrollToTop(500);
	clearForm();
	document.getElementById("rr_raw").value = rrurl;
	if(data == "metapsicov"){
		document.getElementById("srmin").value = 5;
	}
	document.getElementById("contacttypelong").checked = true;
	document.getElementById('submitbutton').focus();
}

function native_analysis(data, prompt, nativeurl){
	var r = confirm("Load data for ".concat(prompt, "?"));
	if (r == false) {
		return;
	}
	scrollToTop(500);
	clearForm();
	document.getElementById("atomtypecb").checked     = true;
	document.getElementById("pdb_raw").value          = nativeurl;
	document.getElementById("native_rr_flag").checked = true;
	document.getElementById('distthresholdrow').style.display = '';
	document.getElementById("rr_raw").value            = "";
	document.getElementById("rrfile").value            = "";
	document.getElementById('tr_rr_raw').style.display = 'none';
	document.getElementById('tr_rr_upload').style.display = 'none';
	if(data == "metapsicov"){
		document.getElementById("srmin").value = 5;
		document.getElementById('contacttypeall').checked = true;
	}
	if(data == "casp"){
		document.getElementById("srmin").value = 6;
		document.getElementById("contacttypelong").checked = true;
	}
	//event.preventDefault();
	document.getElementById('submitbutton').focus();
}

function validateForm(){
	var rrraw = document.getElementById("rr_raw").value;
	var rrfile = document.getElementById("rrfile").value;
	var pdbraw = document.getElementById("pdb_raw").value;
	var pdbfile = document.getElementById("pdbfile").value;
	var na_flag = document.getElementById("native_rr_flag").checked;
	var rrrawstatus = -1;
	var pdrawstatus = -1;
	var rrfilstatus = -1;
	var pdfilstatus = -1;
	if (rrraw == null || rrraw == "") rrrawstatus = 0;
	else rrrawstatus = 1;
	if (rrfile == null || rrfile == "") rrfilstatus = 0;
	else rrfilstatus = 1;
	if (pdbraw == null || pdbraw == "") pdrawstatus = 0;
	else pdrawstatus = 1;
	if (pdbfile == null || pdbfile == "") pdfilstatus = 0;
	else pdfilstatus = 1;
	// at least one input must be there
	if(rrrawstatus == 0 && rrfilstatus == 0 && pdrawstatus == 0 && pdfilstatus == 0){
		alert("Please provide an RR-File or a PDB file!");
		return false;
	}
	// either RR file or RR text field, not both
	if(rrrawstatus == 1 && rrfilstatus == 1){
		alert("Please supply RR file either through the text field or upload a local file! Not Both!");
		return false;
	}
	// either PDB file or PDB text field, not both
	if(pdrawstatus == 1 && pdfilstatus == 1){
		alert("Please supply PDB file either through the text field or upload a local file! Not Both!");
		return false;
	}
	// When RR is supplied, native flag should not be selected
	if((rrrawstatus == 1 || rrfilstatus == 1) && na_flag){
		alert("Please empty the RR file input for PDB analysis!");
		return false;
	}
	// When no RR and pdb is supplied, native flag should be selected
	if((rrrawstatus == 0 && rrfilstatus == 0) && (pdrawstatus == 1 || pdfilstatus == 1) && !na_flag){
		alert("Please check the checkbox for 'Analyze PDB's contacts'! Otherwise, supply RR File for evaluating RR against the native PDB.");
		return false;
	}
	// Values for short-, medium-, and long-range
	var srmin = parseInt(document.getElementById("srmin").value);
	var mrmin = parseInt(document.getElementById("mrmin").value);
	var lrmin = parseInt(document.getElementById("lrmin").value);
	var srmax = parseInt(document.getElementById("srmax").value);
	var mrmax = parseInt(document.getElementById("mrmax").value);
	var lrmax = parseInt(document.getElementById("lrmax").value);
	if(srmin >= srmax){
		alert("Minimum must be less than maximum in sequence separation thresholds for short-range type!");
		return false;
	}
	if(mrmin >= mrmax){
		alert("Minimum must be less than maximum in sequence separation thresholds for medium-range type!");
		return false;
	}
	if(lrmin >= lrmax){
		alert("Minimum must be less than maximum in sequence separation thresholds for long-range type!");
		return false;
	}
	if(srmin < 0){
		alert("Minimum sequence separation threshold for short-range type must be positive!");
		return false;
	}
	if(mrmin - srmax != 1){
		alert("Maximum sequence separation threshold for short-range type must be 1 less than minimum sequence separation threshold for medium-range type!");
		return false;
	}
	if(lrmin - mrmax != 1){
		alert("Maximum sequence separation threshold for medium-range type must be 1 less than minimum sequence separation threshold for long-range type!");
		return false;
	}
}

function scrollToTop(scrollDuration) {
    var scrollStep = -window.scrollY / (scrollDuration / 15),
        scrollInterval = setInterval(function(){
        if ( window.scrollY != 0 ) {
            window.scrollBy( 0, scrollStep );
        }
        else clearInterval(scrollInterval); 
    },15);
}

function toggle_it(itemID){ 
	// Toggle visibility between none and '' 
	if ((document.getElementById(itemID).style.display == 'none')) { 
		document.getElementById(itemID).style.display = '';
		//event.preventDefault();
	}else{ 
		document.getElementById(itemID).style.display = 'none'; 
		//event.preventDefault();
	}    
} 

function check_and_show_dthres_row(){
	if(document.getElementById('native_rr_flag').checked){
		document.getElementById('distthresholdrow').style.display = ''; 
		document.getElementById('tr_rr_raw').style.display = 'none';
		document.getElementById('tr_rr_upload').style.display = 'none';
		document.getElementById("rr_raw").value            = "";
		document.getElementById("rrfile").value            = "";
		//event.preventDefault();
	}else{ 
		document.getElementById('distthresholdrow').style.display = 'none'; 
		document.getElementById('tr_rr_raw').style.display = '';
		document.getElementById('tr_rr_upload').style.display = '';
		//event.preventDefault();
	}
}

function updateProgress(newValue) {
	document.getElementById("progressbar").value = newValue;
}