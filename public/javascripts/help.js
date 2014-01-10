function showValue(newValue, span_id) {
	document.getElementById(span_id).innerHTML=newValue;
}

var addError = function(selector, message){
    $(selector).parents('.form-group').addClass('has-error');
    $(selector).siblings('.help-block').text(message);
}

var removeError = function(selector){
    $(selector).parents('.form-group').removeClass('has-error');
}

$('#submitForm').click(function(event){
    event.preventDefault();
    removeError('#fasta_data');
    removeError('#mut_data');
    removeError('#ps_data');
    valid = true;
    if(!fastaValid()){
        addError('#fasta_data', 'Invalid FASTA format! A more detailed error message to come!');
        valid = false;
    }
    if(!mutValid()){
        addError('#mut_data', 'Invalid mutation format! A more detailed error message to come!');
        valid = false;
    }
    if($('#incl_ps-2').is(':checked') && !psValid()){
        addError('#ps_data', 'Invalid phosphorylation format! A more detailed error message to come!');
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
            