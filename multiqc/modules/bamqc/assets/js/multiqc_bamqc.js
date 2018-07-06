// Set up listeners etc on page load
$(function () {

    // Add the pass / fails counts to each of the BamQC submodule headings
    $.each(bamqc_passfails, function(k, vals){
        var pid = '#bamqc_'+k;
        var total = 0;
        var v = { 'pass': 0, 'fail': 0 };
        var add_warning = false;
        $.each(vals, function(s_name, status){
            if (s_name != 'warning') {
                total += 1;
                v[status] += 1;
            } else if (s_name == 'warning' && status) {
                add_warning = true;
            }
        });
        var pass_percent = (v['pass']/total)*100;
        var fail_percent = (v['fail']/total)*100;
        var pass_number_of_digits = v['pass'].toString().length;
        var fail_number_of_digits = v['fail'].toString().length;
        var pass = pass_percent; // variable used to specify the width of a progress-bar-pass 
        var fail = fail_percent; // vairable used to specify the width of a progress-bar-fail
        // Each digit needs around 8px of space
        // Because the width for either samples that are passed or samples that are failed will be their corresponding percentage times the total width,
        // the total width of a progress bar is 100px, then in order to make all digits visible, a number needs at least (8 * the number of digits it has)% * 100px of space
        if (pass_percent != 0 && pass_percent < pass_number_of_digits * 8) {
            pass = pass_number_of_digits * 8;
            fail = 100 - pass_number_of_digits * 8;
        } else if (fail_percent != 0 && fail_percent < fail_number_of_digits * 8) {
            fail = fail_number_of_digits * 8;
            pass = 100 - fail_number_of_digits * 8;
        }
        var p_bar = '<div class="progress bamqc_passfail_progress"> \
            <div class="progress-bar progress-bar-info" style="width: '+pass+'%" title="'+v['pass']+'&nbsp;/&nbsp;'+total+' samples passed">'+v['pass']+'</div> \
            <div class="progress-bar progress-bar-warning" style="width: '+fail+'%" title="'+v['fail']+'&nbsp;/&nbsp;'+total+' samples are more than two standard deviations away from the mean">'+v['fail']+'</div> \
        </div>';
        $(pid).append(p_bar);
        if (add_warning) {
            var p_warning = '<div class="alert bamqc_alert_warning"> \
            <strong>Warning: </strong><span style="display:inline-block; width:23px;"></span>( mean - 2 standard deviations ) < 0. \
            </div>';
            $(pid).append(p_warning);
        }
    });

    // Create popovers on click
    $('.mqc-section-bamqc .bamqc_passfail_progress .progress-bar').mouseover(function(){
        // Does this element already have a popover?
        if ($(this).attr('data-original-title')) { return false; }
        // Create it
        var pid = $(this).closest('h3').attr('id');
        var k = pid.substr(6);
        var vals = bamqc_passfails[k];
        var passes = $(this).hasClass('progress-bar-info') ? true : false;
        var fails = $(this).hasClass('progress-bar-warning') ? true : false;
        var pclass = '';
        if(passes){ pclass = 'info'; }
        if(fails){ pclass = 'warning'; }
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

    $('html').on('click', function(e) {
        if(!$(e.target).is('.progress-bar') && $(e.target).closest('.popover').length == 0) {
            $('.mqc-section-bamqc .bamqc_passfail_progress .progress-bar').popover('hide');
        }
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
