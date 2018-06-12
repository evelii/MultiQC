// Set up listeners etc on page load
$(function () {

    // Add the pass / fails counts to each of the BamQC submodule headings
    $.each(bamqc_passfails, function(k, vals){
        var pid = '#bamqc_'+k;
        var total = 0;
        var v = { 'pass': 0, 'fail': 0 };
        $.each(vals, function(s_name, status){
            total += 1;
            v[status] += 1;
        });
        var p_bar = '<div class="progress bamqc_passfail_progress"> \
            <div class="progress-bar progress-bar-success" style="width: '+(v['pass']/total)*100+'%" title="'+v['pass']+'&nbsp;/&nbsp;'+total+' samples passed">'+v['pass']+'</div> \
            <div class="progress-bar progress-bar-danger" style="width: '+(v['fail']/total)*100+'%" title="'+v['fail']+'&nbsp;/&nbsp;'+total+' samples are more than two standard deviations away from the mean">'+v['fail']+'</div> \
        </div>';
        $(pid).append(p_bar);
    });

    // Create popovers on click
    $('.mqc-section-bamqc .bamqc_passfail_progress .progress-bar').mouseover(function(){
        // Does this element already have a popover?
        if ($(this).attr('data-original-title')) { return false; }
        // Create it
        var pid = $(this).closest('h3').attr('id');
        var k = pid.substr(6);
        var vals = bamqc_passfails[k];
        var passes = $(this).hasClass('progress-bar-success') ? true : false;
        var fails = $(this).hasClass('progress-bar-danger') ? true : false;
        var pclass = '';
        if(passes){ pclass = 'success'; }
        if(fails){ pclass = 'danger'; }
        var samples = Array();
        $.each(vals, function(s_name, status){
            if(status == 'pass' && passes){ samples.push(s_name); }
            else if(status == 'fail' && fails){ samples.push(s_name); }
        });
        $($(this)).popover({
            animation: true,
            title: $(this).attr('title'),
            content: samples.sort().join('<br>'),
            html: true,
            trigger: 'hover click focus',
            placement: 'bottom auto',
            template: '<div class="popover popover-'+pclass+'" role="tooltip"> \
                <div class="arrow"></div>\
                <h3 class="popover-title"></h3>\
                <div class="bamqc-popover-intro">\
                    Click bar to fix in place <br>\
                    <a href="#" class="bamqc-status-highlight"><span class="glyphicon glyphicon-pushpin"></span> Highlight these samples</a><br>\
                    <a href="#" class="bamqc-status-hideothers"><span class="glyphicon glyphicon-eye-close"></span> Show only these samples</a>\
                </div>\
                <div class="popover-content"></div>\
            </div>'
        }).popover('show');
    });

    // Listener for Status higlight click
    $('.mqc-section-bamqc .bamqc_passfail_progress').on('click', '.bamqc-status-highlight', function(e){
        e.preventDefault();
        // Get sample names and highlight colour
        var samples = $(this).parent().parent().find('.popover-content').html().split('<br>');
        var f_col = mqc_colours[mqc_colours_idx];
        // Add sample names to the toolbox
        for (i = 0; i < samples.length; i++) {
            var f_text = samples[i];
            $('#mqc_col_filters').append('<li style="color:'+f_col+';"><span class="hc_handle"><span></span><span></span></span><input class="f_text" value="'+f_text+'"/><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>');
        }
        // Apply highlights and open toolbox
        apply_mqc_highlights();
        mqc_toolbox_openclose('#mqc_cols', true);
        // Update next highlight colour
        mqc_colours_idx += 1;
        if(mqc_colours_idx >= mqc_colours.length){ mqc_colours_idx = 0; }
        $('#mqc_colour_filter_color').val(mqc_colours[mqc_colours_idx]);
        // Hide the popover
        $(this).closest('.popover').popover('hide');
    });

    // Listener for Status hide others click
    $('.mqc-section-bamqc .bamqc_passfail_progress').on('click', '.bamqc-status-hideothers', function(e){
        e.preventDefault();
        // Get sample names
        var samples = $(this).parent().parent().find('.popover-content').html().split('<br>');
        // Check if we're already hiding anything, remove after confirm if so
        if($('#mqc_hidesamples_filters li').length > 0){
            if(!confirm($('#mqc_hidesamples_filters li').length+' Hide filters already exist - discard?')){
                return false;
            } else {
                $('#mqc_hidesamples_filters').empty();
            }
        }
        // Set to "show only" and disable regex
        $('.mqc_hidesamples_showhide[value="show"]').prop('checked', true);
        $('#mqc_hidesamples .mqc_regex_mode .re_mode').removeClass('on').addClass('off').text('off');
        // Add sample names to the toolbox
        for (i = 0; i < samples.length; i++) {
            var f_text = samples[i];
            $('#mqc_hidesamples_filters').append('<li><input class="f_text" value="'+f_text+'" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>');
        }
        // Apply highlights and open toolbox
        apply_mqc_hidesamples();
        mqc_toolbox_openclose('#mqc_hidesamples', true);
        // Hide the popover
        $(this).closest('.popover').popover('hide');
    });

});
