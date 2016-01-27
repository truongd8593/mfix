var mfixrunning = false;

$(document).ready(function(){
    $("#curlstartstop").text('curl -X PUT '+document.location.origin+'/stop');
    updateCurlCommands();
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
                updateCurlCommands();
            }
        });
    });

    $('#get').click(function(){
        var varname = ["mfix",
                           $("#getmodname").val(),
                           $("#getvarname").val()].join('.');
        var value = $("#getvalue").val();

        $.ajax({
            url: 'get/'+varname,
            type: 'GET',
            success: function(response) {
                $("#getresponse").text(response);
                $("#status").text('Successfully got value');
                $("#status").css('color', 'green');
            },
            error: function(response) {
                $("#getresponse").text(response);
                $("#status").text('Error retrieving value');
                $("#status").css('color', 'red');
            }
        });
    });

    $('#set').click(function(){
        var varname = ["mfix",
                           $("#setmodname").val(),
                           $("#setvarname").val()].join('.');
        var value = $("#setvalue").val();

        $.ajax({
            url: 'set/'+varname,
            type: 'POST',
            data: {'varvalue':value},
            success: function(response) {
                $("#status").text('Successfully set value');
                $("#status").css('color', 'green');
            },
            error: function(response) {
                $("#status").text('Error doing set value');
                $("#status").css('color', 'red');
            }
        });
    });

    $('#step').click(function(){

        var stepcount = $("#stepcount").val();
        $.ajax({
            url: 'step',
            type: 'POST',
            data: {'stepcount':stepcount},
            success: function(response) {
                $("#status").text('Successfully did timesteps');
                $("#status").css('color', 'green');
            },
            error: function(response) {
                $("#status").text('Error doing timesteps');
                $("#status").css('color', 'red');
            }
        });
    });

    $('#backupres').click(function(){

        $.ajax({
            url: 'backupres',
            type: 'POST',
            success: function(response) {
                $("#status").text('Successfully backed up resource files');
                $("#status").css('color', 'green');
            },
            error: function(response) {
                $("#status").text('Error while backing up resource files');
                $("#status").css('color', 'red');
            }
        });
    });

    $('#reinit').click(function(){

        $.ajax({
            url: 'reinit',
            type: 'POST',
            success: function(response) {
                $("#status").text('Successfully reinitialized');
                $("#status").css('color', 'green');
            },
            error: function(response) {
                $("#status").text('Reinitialize failed');
                $("#status").css('color', 'red');
            }
        });
    });

    $("input, select").change(updateCurlCommands);
    $("input").keydown(updateCurlCommands);
    $("input").keyup(updateCurlCommands);

    $('.stepcount').keyup(function () {
        if (this.value != this.value.replace(/[^0-9]/g, '')) {
            this.value = this.value.replace(/[^0-9]/g, '');
        }
    });

});

function updateCurlCommands() {
    $("#curlbackupres").text('curl -X POST '+document.location.origin+'/backupres');
    $("#curlreinit").text('curl -X POST '+document.location.origin+'/reinit');
    $("#curlstep").text('curl -X PUT '+document.location.origin+'/step');

    var varname = ["mfix",
                       $("#setmodname").val(),
                       $("#setvarname").val()].join('.');
    var value = $("#setvalue").val();
    $("#curlset").text('curl -X POST '+document.location.origin+'/set/'+varname+' -d "varvalue='+value+'"');

    var varname = ["mfix",
                       $("#getmodname").val(),
                       $("#getvarname").val()].join('.');
    $("#curlget").text('curl -X GET '+document.location.origin+'/get/'+varname);

    if (mfixrunning) {
        $("#running").text('MFIX IS RUNNING');
        $("#curlstartstop").text('curl -X PUT '+document.location.origin+'/stop');
        $('#step').prop('disabled',true);
    } else {
        $("#running").text('MFIX IS STOPPED');
        $("#curlstartstop").text('curl -X PUT '+document.location.origin+'/start');
        $('#step').prop('disabled',false);
    }

}

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
