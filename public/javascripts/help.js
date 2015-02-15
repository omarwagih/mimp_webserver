Array.prototype.unique = function() {
    var o = {}, i, l = this.length, r = [];
    for(i=0; i<l;i+=1) o[this[i]] = this[i];
    for(i in o) r.push(o[i]);
    return r;
};

function showValue(newValue, span_id) {
	document.getElementById(span_id).innerHTML=newValue;
}

function toggleAdvanced(){
    if($("#advanced-div").is(":visible")){

        $("#adv-label").text("Show advanced options");
        $("#adv-icon").removeClass("fa-chevron-down");
        $("#adv-icon").addClass("fa-chevron-right");
        $("#advanced-div").slideUp();
    }else{
        $("#adv-label").text("Hide advanced options");
        $("#adv-icon").removeClass("fa-chevron-right");
        $("#adv-icon").addClass("fa-chevron-down");
        $("#advanced-div").slideDown();
    }

}

var addError = function(selector, message){
    $(selector).parents('.form-group').addClass('has-error');
    $(selector).siblings('.help-block').html(message);
}

var addWarning = function(selector, message){
    $(selector).parents('.form-group').addClass('has-warning');
    $(selector).siblings('.help-block').html(message);
}

var removeError = function(selector){
	$(selector).parents('.form-group').find('.help-block').empty();
    $(selector).parents('.form-group').removeClass('has-error').removeClass('has-warning');;
}

$('#model_data').change(function(){
	v = $( "#model_data option:selected" ).val();
	if(v == "lconf"){
		addWarning("#model_data", "Predicted models are purely speculative and should be used with caution. They are based on the kinase-substrate relationships and not kinase-peptide. We strongly urge users who use this dataset to manually verify that the logos appear valid before proceeding with downstream analyses. For more information on how these specifcity models are predicted, please read the <a class='warning-link' href='../help'>help section</a>.")
	}else{
		removeError("#model_data");
	}
});

$('#submitForm').click(function(event){
    event.preventDefault();
    removeError('#fasta_data');
    removeError('#mut_data');
    removeError('#ps_data');
    valid = true;
    
    mut_data = $('#mut_data').val().split('\n').unique();
    MUT_LIMIT = 2000;
    if(mut_data.length > MUT_LIMIT){
    	jobidfield.show();
    	jobidfield.find('label').css('visibility', 'hidden');
    	jobidfield.find('.help-block').html('You have provided over '+MUT_LIMIT.toLocaleString()+' mutations. The web interface was designed to handle smaller analyses. To process a large number of mutations, please download the MIMP R package from the <a target="_blank" class="error-link" href="/download">download page</a> and use it locally.')
    	jobidfield.addClass('has-error');
    	window.scrollTo(0, 0);
    	return;
    }
    // Set mutation data to be the unique ones
    $('#mut_data').val(mut_data.join('\n'));
    
    fv = fastaValid();
    if(!fv.valid){
        addError('#fasta_data', 'Error on line '+fv.line+': '+fv.msg+'!<br>Please provide sequences in correct FASTA format. For more information see the <a class="error-link" target="_blank" href="/help">help page</a> or click "Load sample" for example data.');
        valid = false;
    }
    if(!mutValid()){
        addError('#mut_data', 'Invalid mutation data! Please provide mutation data, line separated. For example, TP53 R282W. For more information see the <a class="error-link" target="_blank" href="/help">help page</a> or click "Load sample" for example data.');
        valid = false;
    }
    if($('#incl_ps-2').is(':checked') && !psValid()){
        addError('#ps_data', 'Invalid phosphosite data! Please provide phosphosite data, line separated. For example, TP53 150. For more information see the <a class="error-link" target="_blank" href="/help">help page</a>.');
        valid = false;
    }
    if(!valid) return;
    
    $('form').submit();
});

var fastaValid = function(){
    s = $('#fasta_data').val();
    var aa_regex = new RegExp('[^A-Z\*]+','gi');
    var header_regex = new RegExp('^>[\\w].*');
    s = s.trim().split('\n');
    header = false;
    seq = 0;
    if(s[0][0] != '>'){
    	return {valid:false, msg:"headers must begin with '>'", line:1};
    	//return false;
    }
    for(var i = 0; i < s.length; i++){
        z = s[i];
        if(!z.trim()){
            //console.log(i + ' blank');
            continue;
        }else if(z[0] == '>'){ 
        	if(!header_regex.test(z)){
            	return {valid:false, msg:"headers must begin with '>', followed by alphanumerics", line:i+1};
        		// Header should start with a letter (no space);
        		// return false;
        	}
            if( i != 0 & seq == 0 ){
            	return {valid:false, msg:"header has no sequence", line:i+1};
                //console.log(i + ' some header has no sequences');
                //return false;
            }
            seq = 0;
            //console.log(i + ' header');
            header = true;
        }else if(aa_regex.test(z)){
            //console.log(i + ' invalid line');
            //Invalid fasta
        	return {valid:false, msg:"invalid amino acid sequence", line:i+1};
            //return false;
        }else{
            //console.log(i + ' incr');
            seq++;
        }
    }
    //console.log(seq)
    if(!header){
    	return {valid:false, msg:"sequence has no header", line:s.length};
    }
    if(seq < 1){
    	return {valid:false, msg:"header has no sequence", line:s.length};
        //return false;
    }

	return {valid:true, msg:"", line:0};
    //return true;

}

            
function mutValid(){
    var reg = /^\S+\s[ACDEFGHIKLMNPQRSTVWY]\d+[ACDEFGHIKLMNPQRSTVWY]$/i
    var z = $('#mut_data').val()

    if(!z.trim()) return false;

    z = z.split("\n");
    for(var i=0; i < z.length; i++){
        console.log(z[i]);
        if(!z[i]){
            continue;
        }else if(!reg.test(z[i])){
            return false;
        }
    }
    return true;
}

function psValid(){
    var reg = /^\S+\s\d+$/i
    var z = $('#ps_data').val()

    if(!z.trim()) return false;
    z = z.split("\n");
    for(var i=0; i < z.length; i++){
        console.log(z[i]);
        if(!z[i]){
            continue;
        }if(!reg.test(z[i])){
            return false;
        }
    }
    return true;
}
            