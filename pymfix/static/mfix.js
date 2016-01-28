var mfixrunning = false;
var req_common = '>>> import requests; requests.post("'+document.location.origin;

$(document).ready(function(){
    $("#curlstartstop").text(req_common+'/stop")');
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

    $('#pausetime').change(function(){
        var pausetime = $("#pausetime").val();
        $.ajax({
            url: 'set/mfix.main.pausetime',
            type: 'POST',
            data: {'varvalue':pausetime},
            success: function(response) {
                $("#status").text('Successfully set pausetime');
                $("#status").css('color', 'green');
            },
            error: function(response) {
                $("#status").text('Error setting pausetime');
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

    $('#exit').click(function(){
        $.ajax({
            url: 'exit',
            type: 'POST',
            success: function(response) {
                $("#status").text('Successfully exited');
                $("#status").css('color', 'green');
            },
            error: function(response) {
                $("#status").text('Exit failed');
                $("#status").css('color', 'red');
            }
        });
    });

    $('#abort').click(function(){
        $.ajax({
            url: 'abort',
            type: 'POST',
            success: function(response) {
                $("#status").text('Successfully aborted');
                $("#status").css('color', 'green');
            },
            error: function(response) {
                $("#status").text('Abort failed');
                $("#status").css('color', 'red');
            }
        });
    });

    $("input, select").change(updateCurlCommands);
    $("input").keydown(updateCurlCommands);
    $("input").keyup(updateCurlCommands);

    // TODO: better integer validation
    $('#stepcount').keyup(function () {
        if (this.value != this.value.replace(/[^0-9]/g, '')) {
            this.value = this.value.replace(/[^0-9]/g, '');
        }
    });

    // TODO: better float validation
    $('#pausetime').keyup(function () {
        if (this.value != this.value.replace(/[^0-9\.]/g, '')) {
            this.value = this.value.replace(/[^0-9\.]/g, '');
        }
    });

});

function updateCurlCommands() {
    $("#curlbackupres").text(req_common+'/backupres")');
    $("#curlreinit").text(req_common+'/reinit")');
    $("#curlexit").text(req_common+'/exit")');
    $("#curlabort").text(req_common+'/abort")');
    $("#curlpausetime").text(req_common+'/set/mfix.main.pausetime", data={"varvalue":"'+$("#pausetime").val()+'"})');
    $("#curlstep").text(req_common+'/step'+'", data={"stepcount":"'+$("#stepcount").val()+'"})');

    var varname = ["mfix",
                   $("#setmodname").val(),
                   $("#setvarname").val()].join('.');
    var value = $("#setvalue").val();
    $("#curlset").text(req_common+'/set/'+varname+'", data={"varvalue":"'+value+'"})');

    var varname = ["mfix",
                   $("#getmodname").val(),
                   $("#getvarname").val()].join('.');
    $("#curlget").text(req_common+'/get/'+varname+'").text');

    if (mfixrunning) {
        $("#running").text('MFIX IS RUNNING');
        $("#curlstartstop").text(req_common+'/stop")');
        $('#step').prop('disabled',true);
        $('#pausetime').prop('disabled',true);
    } else {
        $("#running").text('MFIX IS STOPPED');
        $("#curlstartstop").text(req_common+'/start")');
        $('#step').prop('disabled',false);
        $('#pausetime').prop('disabled',false);
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
