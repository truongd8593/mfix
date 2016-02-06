var mfixrunning = false;
var req_common = '>>> import requests; requests.post("'+document.location.origin;

$(document).ready(function(){

    $(".notice").hide();
    $('.notice').click(function(){
        $(this).fadeOut();
    });

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
                $(".notice").hide();
                $("div.success").text('Successfully got value');
                $("div.success").fadeIn();
            },
            error: function(response) {
                $("#getresponse").text("");
                $(".notice").hide();
                $("div.error").text('Error retrieving value');
                $("div.error").fadeIn();
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
                $(".notice").hide();
                $("div.success").text('Successfully set value');
                $("div.success").fadeIn();
            },
            error: function(response) {
                $(".notice").hide();
                $("div.error").text('Error setting value');
                $("div.error").fadeIn();
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
                $(".notice").hide();
                $("div.success").text('Successfully did timestep(s)');
                $("div.success").fadeIn();
            },
            error: function(response) {
                $(".notice").hide();
                $("div.error").text('Error doing timestep(s)');
                $("div.error").fadeIn();
            }
        });
    });

    $('#write_dbg_vt').click(function(){
        $.ajax({
            url: 'write_dbg_vt',
            type: 'POST',
            success: function(response) {
                $(".notice").hide();
                $("div.success").text('Successfully called WRITE_DBG_VTU_AND_VTP_FILES');
                $("div.success").fadeIn();
            },
            error: function(response) {
                $(".notice").hide();
                $("div.error").text('Error while calling WRITE_DBG_VTU_AND_VTP_FILES');
                $("div.error").fadeIn();
            }
        });
    });

    $('#backupres').click(function(){
        $.ajax({
            url: 'backupres',
            type: 'POST',
            success: function(response) {
                $(".notice").hide();
                $("div.success").text('Successfully backed up resource files');
                $("div.success").fadeIn();
            },
            error: function(response) {
                $(".notice").hide();
                $("div.error").text('Error while backing up resource files');
                $("div.error").fadeIn();
            }
        });
    });

    $('#reinit').click(function(){
        $.ajax({
            url: 'reinit',
            type: 'POST',
            success: function(response) {
                $(".notice").hide();
                $("div.success").text('Successfully reinitialized');
                $("div.success").fadeIn();
            },
            error: function(response) {
                $(".notice").hide();
                $("div.error").text('Reinitialize failed');
                $("div.error").fadeIn();
            }
        });
    });

    $('#exit').click(function(){
        $.ajax({
            url: 'exit',
            type: 'POST',
            success: function(response) {
                $(".notice").hide();
                $("div.success").text('Successfully exited');
                $("div.success").fadeIn();
            },
            error: function(response) {
                $(".notice").hide();
                $("div.error").text('Exit failed');
                $("div.error").fadeIn();
            }
        });
    });

    $('#abort').click(function(){
        $.ajax({
            url: 'abort',
            type: 'POST',
            success: function(response) {
                $(".notice").hide();
                $("div.error").text('Successfully aborted');
                $("div.error").fadeIn();
            },
            error: function(response) {
                $(".notice").hide();
                $("div.error").text('Abort failed');
                $("div.error").fadeIn();
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
});

function updateCurlCommands() {
    $("#curlwritedbgvt").text(req_common+'/write_dbg_vt")');
    $("#curlbackupres").text(req_common+'/backupres")');
    $("#curlreinit").text(req_common+'/reinit")');
    $("#curlexit").text(req_common+'/exit")');
    $("#curlabort").text(req_common+'/abort")');
    $("#curlstep").text(req_common+'/step'+'", data={"stepcount":"'+$("#stepcount").val()+'"})');

    var varname = ["mfix",
                   $("#setmodname").val(),
                   $("#setvarname").val()].join('.');
    var value = $("#setvalue").val();
    $("#curlset").text(req_common+'/set/'+varname+'", data={"varvalue":"'+value+'"})');

    var varname = ["mfix",
                   $("#getmodname").val(),
                   $("#getvarname").val()].join('.');
    $("#curlget").text(req_common.replace("requests.post","requests.get")+'/get/'+varname+'").text');

    if (mfixrunning) {
        $("#running").text('MFIX IS RUNNING');
        $("#curlstartstop").text(req_common+'/stop")');
        $('button, input, select').prop('disabled',true);
    } else {
        $("#running").text('MFIX IS STOPPED');
        $("#curlstartstop").text(req_common+'/start")');
        $('button, input, select').prop('disabled',false);
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
