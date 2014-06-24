Array.prototype.unique = function() {
    var o = {}, i, l = this.length, r = [];
    for(i=0; i<l;i+=1) o[this[i]] = this[i];
    for(i in o) r.push(o[i]);
    return r;
};

function showValue(newValue, span_id) {
	document.getElementById(span_id).innerHTML=newValue;
}

var addError = function(selector, message){
    $(selector).parents('.form-group').addClass('has-error');
    $(selector).siblings('.help-block').html(message);
}

var removeError = function(selector){
	$(selector).parents('.form-group').find('.help-block').empty();
    $(selector).parents('.form-group').removeClass('has-error');
}

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
    
    if(!fastaValid()){
        addError('#fasta_data', 'Invalid data! Please provide sequences in FASTA format. For more information see the <a class="error-link" target="_blank" href="/help">help page</a> or click "Load sample" for example data.');
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
    var aa_regex = new RegExp('[^ARNDCQEGHILKMFPSTWYVBZ\*_-]+','gi');
    s = s.trim().split('\n');
    header = false;
    seq = 0;
    if(s[0][0] != '>'){
        return false;
    }
    for(var i = 0; i < s.length; i++){
        z = s[i];
        if(!z.trim()){
            //console.log(i + ' blank');
            continue;
        }else if(z[0] == '>'){ 
            if( i != 0 & seq == 0 ){
                //console.log(i + ' some header has no sequences');
                return false;
            }
            seq = 0;
            //console.log(i + ' header');
            header = true;
        }else if(aa_regex.test(z)){
            //console.log(i + ' invalid line');
            //Invalid fasta
            return false;
        }else{
            //console.log(i + ' incr');
            seq++;
        }
    }
    //console.log(seq)
    if(seq < 1 || !header){
        return false;
    }
    return true;

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
            