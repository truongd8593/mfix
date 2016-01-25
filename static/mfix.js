$(document).ready(function(){
    var mfixrunning = true;
    $('a.toggler').click(function(){

        if (mfixrunning) {
            url = "stop"
        } else {
            url = "start"
        }
        $(this).toggleClass('off');
        $.ajax({
            url: url,
            type: 'PUT',
            success: function(response) {
                mfixrunning = !mfixrunning;
                if (mfixrunning) {
                    $("#running").text('MFIX IS RUNNING');
                    $("#selectable").text('curl -X PUT http://localhost:5000/stop');
                } else {
                    $("#running").text('MFIX IS STOPPED');
                    $("#selectable").text('curl -X PUT http://localhost:5000/start');
                }
            }
        });
    });

    $('#set').click(function(){

        $.ajax({
            url: 'set',
            type: 'POST',
            success: function(response) {
                var retVal = true
                if (retVal) {
                    $("#status").text('Successfully set value');
                    $("#status").css('color', 'green');
                } else {
                    $("#status").text('Error while setting value');
                    $("#status").css('color', 'red');
                }
            }
        });
    });
});

function selectText(elem) {
    if (document.selection) {
        var range = document.body.createTextRange();
        range.moveToElementText(elem);
        range.select();
    } else if (window.getSelection) {
        var range = document.createRange();
        range.selectNode(elem);
        window.getSelection().addRange(range);
    }
}
