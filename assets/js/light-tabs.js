(function($) {

	var tabs =  $(".tabs li a");
  
	tabs.click(function() {
		var content = this.hash.replace('/','');
		tabs.removeClass("active");
		$(this).addClass("active");
    $("#tab-content").find('p').hide();
    $(content).fadeIn(200);
	});

})(jQuery);
