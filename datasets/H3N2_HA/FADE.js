var _colorizerB = d3.interpolateRgb  (d3.rgb(0,0,255),d3.rgb(255,255,255));
var _colorizerR = d3.interpolateRgb  (d3.rgb(255,255,255),d3.rgb(255,0,0));
var _use_q_values = false;
var selected = true;

function dict_to_array (dict) {
        ar = []
        for (k in dict) {
            ar.push (dict[k])
        }
        return ar;
    }
    
//For displaying table with Posteriors
function display_column_map (row) { 
    //result = [parseInt(row[0]),row[3]];
    result = [parseInt(row[0])];
	
	for (k = 1; k < 21; k++) {
        result.push((row[k]));
    }
    return result;
}

//For displaying table with BFs
function display_column_map_q (row) { 
     //result = [parseInt(row[0]),row[3]];
    result = [parseInt(row[0])];
	
	for (k = 1; k < 21; k++) {
        result.push((row[k]));
    }
    return result;
}

//For setting up plot data
function property_for_plotting (row) { 
    result = [row[0]];
    for (k = 1; k <= 20; k++) {
        result.push(row[k]);
    }
    return result;
}

function row_display_filter (d) {

//Any row, with at least one val > thres must get displayed. Any elements greater must be in red. 
   // if (d.slice(2).reduce (function (a,b) {return a+b;}) == 0.0) {return false;} 
    //console.log (d, this);
    for (k=1;k<21;k++) {
        if (d[k] > this) return true;
    } 
    return false;
};

function initial_display  () {
    load_data_summary ("summary_div", data_summary);
    $('#filter_on_pvalue').trigger ('submit');
	
    $('#property_plot_collapse').on('show', function () {
        plot_property_graphs("property_plot_svg",FADE_posteriors); //Using a matrix from html
    });
}

function set_handlers  (file_id) {
    $('body').attr ('data-job-id', file_id);
    $('#filter_on_pvalue').submit(function (e) {
            cutoff = parseFloat($('#pvalue')[0].value);
            if (_use_q_values) {
                found = load_analysis_results('prime_table',prime_headers, FADE_BFs,display_column_map_q,row_display_filter);            
            } else {
                found = load_analysis_results('prime_table',prime_headers, FADE_posteriors,display_column_map,row_display_filter);
            }
            d3.select ("#total_sites_found").selectAll("span").data (found).html (function (d) {return d;});
            return false;
        });
	
		
	
    $('#site_rate_display').on ('show', (function (e) {
            //alert ("Show");
            //console.log (this);
            return true;
        }));

    
      $( 'body' ).on( 'click', '[data-toggle="modal"]', function(event) {
            display_site_properties ($(this).attr ("data-codon-id"));
    });
 

   $( '#set-p-value' ).click(function(event) {
         d3.select ("#pq_selector").html ("Posterior <span class='caret'></span>");
         _use_q_values = false;
         event.preventDefault(); 
   } );    

   $( '#set-q-value' ).click(function(event) {
         d3.select ("#pq_selector").html ("BF <span class='caret'></span>");
         _use_q_values = true;
         event.preventDefault(); 
   } );    
         
   $( 'body' ).on( 'click', '#property_selector .btn', function(event) {
        event.stopPropagation(); // prevent default bootstrap behavior
        if( $(this).attr('data-toggle') != 'button' ) { // don't toggle if data-toggle="button"
            $(this).toggleClass('active');
        }
         toggle_view("property_plot_svg", parseInt($(this).attr( 'data-property-id' )), $(this).hasClass( 'active' ) ); // button state AFTER the click
    });
	
	$( '#deselectall' ).click(function(event) {
        
		
		
		for (var i=1;i<=20;i++)
		{ 
		  if (selected)
		  {
			 toggle_view("property_plot_svg", i, false ); // button state AFTER the click
			 $("#deselectall").html("Select All");
			 $("#show_property"+i).removeClass("active");
			
		  }
		  else
		  {
		   toggle_view("property_plot_svg", i, true ); // button state AFTER the click
		 $("#deselectall").html("Deselect All");
		 $("#show_property"+i).addClass("active");

		  }
		 
		 }
		 
		 selected = !(selected);
		 
   } );   
   
}

property_plot_done = false;


function display_site_properties  (codon_id) {
    job_id = $('body').attr ('data-job-id');
    url = "/cgi-bin/datamonkey/wrapHyPhyBF.pl?file=prime_site&mode=1&arguments=" + job_id + "-" + codon_id;
    d3.json (url, function (json) { set_exch_matrix (["ACGILMPSTV", "DENQ", "FWY", "HKR"],json, codon_id) });
}

function toggle_view  (property_plot, group, show_hide) {
    if (show_hide) {
        prop = 'visible';
    } else {
        prop = 'hidden';
    }
    d3.select("#"+property_plot).selectAll(".dot"+group)
                  .style("visibility", prop);
}


function plot_property_graphs  (property_plot, property_info) {
    if (!property_plot_done) {
        property_info = property_info.map (property_for_plotting);
        
        property_plot_done = true;
        site_count = property_info.length;
               
        //console.log (d3.extent (property_info.map(function (d){return d[0];})));
                
        var margin = {top: 20, right: 40, bottom: 30, left: 40},
            width  = 940 - margin.left - margin.right,
            height = 500 - margin.top - margin.bottom;

        var x = d3.scale.linear()
            .range([0, width]);

        var y = d3.scale.linear()
            .range([height, 0]);

        var color = d3.scale.category20();

        var xAxis = d3.svg.axis()
            .scale(x)
            .orient("bottom");

        var yAxis = d3.svg.axis()
            .scale(y)
            .orient("left");
 
         var yAxis2 = d3.svg.axis()
            .scale(y)
            .orient("right");
           
        function make_x_axis() {        
            return d3.svg.axis()
                .scale(x)
                 .orient("bottom")
                 .ticks(20);
        }

        function make_y_axis() {        
            return d3.svg.axis()
                .scale(y)
                .orient("left")
                .ticks(20);
        }

        var svg = d3.select("#"+property_plot)
            .attr("width", width + margin.left + margin.right)
            .attr("height", height + margin.top + margin.bottom)
             .append("g")
            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
                
        x.domain ([1,site_count]);
        y.domain ([0,1]);

        svg.append("g")
              .attr("class", "x axis")
              .attr("transform", "translate(0," + height + ")")
              .call(xAxis)
            .append("text")
              .attr("class", "label")
              .attr("x", width)
              .attr("y", 30)
              .style("text-anchor", "end")
              .text("Site index");

        svg.append("g")         
                .attr("class", "grid")
                .call(make_y_axis()
                    .tickSize(-width, 0, 0)
                    .tickFormat("")
                );
        
        svg.append("g")         
                .attr("class", "grid")
                .attr("transform", "translate(0," + height + ")")
                .call(make_x_axis()
                    .tickSize(-height, 0, 0)
                    .tickFormat("")
                );
        
        svg.append("g")
              .attr("class", "y axis")
              .call(yAxis)
            .append("text")
              .attr("class", "label")
              .attr("transform", "rotate(-90)")
              .attr("y", -37)
              .attr("dy", ".71em")
              .style("text-anchor", "end")
              .text("P(Bias>1)");  

        var y2= svg.append("g")
              .attr("class", "y axis")
              .attr("transform", "translate("+width+",0)")
              .call(yAxis2.tickFormat (""));
              
        y2.append("text")
              .attr("class", "label")
              .attr("transform", "rotate(-90)")
              .attr("y", 10)
              .attr("dy", ".71em")
              .style("text-anchor", "end")
              .text("High Posteriors");  

        y2.append("text")
              .attr("class", "label")
              .attr("transform", "rotate(-90)")
              .attr("y", 10)
              .attr("x", -height)
              .attr("dy", ".71em")
              .style("text-anchor", "start")
              .text("Low Posteriors");  

        var legend = svg.selectAll(".legend")
              .data(color.domain())
            .enter().append("g")
              .attr("class", "legend")
              .attr("transform", function(d, i) { return "translate(0," + i * 20 + ")"; });
		
		
		var h = new Object(); //Hash of numbers -> AA names for labels
		h[1] = "Alanine";
		h[2] = "Cysteine";
		h[3] = "Aspartic acid";
		h[4] = "Glutamic acid";
		h[5] = "Phenylalanine";
		h[6] = "Glycine";
		h[7] = "Histidine";
		h[8] = "Isoleucine";
		h[9] = "Lysine";
		h[10] = "Leucine";
		h[11] = "Methionine";
		h[12] = "Asparagine";
		h[13] = "Proline";
		h[14] = "Glutamine";
		h[15] = "Arginine";
		h[16] = "Serine";
		h[17] = "Threonine";
		h[18] = "Valine";
		h[19] = "Tryptophan";
		h[20] = "Tyrosine";
		
		

        for (series = 1; series <= 20; series++) {
            svg.selectAll(".dot"+series)
                  .data(property_info)
                .enter().append("circle")
                  .attr("class", "dot"+series)
                  .attr("r", function (d) {if (d[series] == 0) return 1; return 3.5;})
                  .attr("cx", function(d) { return x(d[0]); })
                  .attr("cy", function(d) { return y(d[series]); })
                  .style("fill", function(d) { return color(series); })
				  .style("opacity", 0.5)
                  .append ("title").text (function (d) { return "Site " + d[0] + ", " + h[series] + " P(Beta>1) =" + d[series];});
            d3.select ("#show_property" + series).style ("color", function(d) { return color(series); }); //Colour buttons on HTML
        }
      
     } 
        
}

function load_data_summary  (id, json_object) {
    //if (error) return console.warn(error);
    reportable_entries = []
    warnings = []
    for (k in json_object) {
        if (k == 'Warnings') {
            warnings = json_object[k];
        } else {
            reportable_entries.push ({"Phase": k, "Information": json_object[k]});
        }
    }
    info_container = d3.select ("#" + id);
    items = info_container.selectAll("dt").data(reportable_entries);
    items.enter().append("dt").text (function (d) {return d["Phase"]});
    item_count = items[0].length;
    current_child = 2;
    for (z = 1; z <= item_count; z++) { 
        info = dict_to_array(reportable_entries[z-1]["Information"]);
        for (y = 0; y < info.length; y++) {
            info_container.insert ("dd",":nth-child("+current_child+")").data([info[y]]).text (function (d) {return d;});
        }
        key = reportable_entries[z-1]["Phase"];
        if (key in warnings) {
            current_child++;
            info_container.insert ("dd",":nth-child("+current_child+")").selectAll("div").data([warnings[key]]).enter().append("div").classed("alert",true).html (function (d) {return d;})
        }
        current_child += info.length+1;
    }
    return 0;
}
    
function load_analysis_results  (id, headers, matrix, column_selector, condition) {    
    header_element = d3.select ('#' + id).select("thead");
    header_element.selectAll("th").remove();
    header_element.selectAll("th").data (headers[0]).enter().append("th").html (function (d,i) //Get header of table
        {return "<a href='#' data-toggle='tooltip' data-placement = 'right' data-html = true title data-original-title='" + headers[1][i] + "'>" + d + "</a>";}
    );
        
    parent_element = d3.select ('#' + id).select("tbody");
    parent_element.selectAll("tr").remove();
    filtered_matrix = matrix.map (column_selector).filter(condition,cutoff); //Get the columns to display in table
    rows = parent_element.selectAll("tr").data(function (d) {return filtered_matrix;});
    conserved = 0;
    changing  = 0;
    rows.enter().append ("tr").selectAll("td").data (function (d) {return d;}).enter().append("td").
        html (function (d,i) {
            if (i && _use_q_values == false) {return ("" + (  (parseFloat(d)) )).substring(1,4);} else {return ("" + (  (parseFloat(d)) )).substring(0,4);}
            return "<b>" + d + "</b> <a href='#site_rate_display' data-toggle='modal' data-codon-id = '" + d + "' data-placement = 'bottom'><i class='icon-th'></i></a>";}
            ).
            classed ("btn-danger", function (d,i,j)  {if (d >= cutoff &&  i>=1 ) {conserved++; return true;} return false;}).
						
						
						
            classed ("btn-success", function (d,i,j) {if (i>1 && (i%2==1 && d <= cutoff && filtered_matrix[j][i-1] < 0 || i%2==0 && filtered_matrix[j][i+1] <= cutoff && d < 0)) {changing++; return true;} return false;});
            
    d3.select ('#' + id).classed ("table-striped table-hover", true);
    $('a').tooltip(); 
    return [filtered_matrix.length, conserved, changing/2];
}

function set_exch_matrix (aa_binning, values, codon_id) {    
    
        var amino_acid_ordering = aa_binning.join("");
        var matrix = [];       
        var unique_rates = [];
        
        for (k = 0; k < 20; k += 1) {
            matrix.push ([]);
            for (k2 = 0; k2 < 20; k2+=1) {
                aa_pair = amino_acid_ordering [k] + amino_acid_ordering [k2];
                if (aa_pair in values) {
                    matrix[k].push (values[aa_pair]);
                    unique_rates.push (values[aa_pair]);
                 } else {
                    aa_pair = amino_acid_ordering [k2] + amino_acid_ordering [k];
                    if (aa_pair in values) {
                        matrix[k].push (values[aa_pair]);
                    } else {
                        matrix[k].push ("");
                        unique_rates.push (values[aa_pair]);
                    }
                }
            }
        }
        
        //console.log (d3.extent (unique_rates), d3.mean (unique_rates), d3.median (unique_rates)); 
        
        var margin = {top: 20, right: 20, bottom: 20, left: 20},
            width  = 400 - margin.left - margin.right,
            height = 400 - margin.top - margin.bottom,
            cellSize = 17;
            
        d3.select ("#site_rate_display_header").html ("Inferred exchangeabilities at codon " + codon_id);
        
        d3.select ("#exch_matrix").selectAll ("g").remove();
        
        var svg = d3.select ("#exch_matrix")
            .attr("width", width + margin.left + margin.right)
            .attr("height", height + margin.top + margin.bottom)
            .append("g")
            .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
                  
      
        var x = d3.scale.ordinal()
            .rangePoints([0, cellSize*20],  1);
            
        x.domain(amino_acid_ordering) ;
      
        var xAxis = d3.svg.axis()
          .scale(x)
          .orient("top");
          
        var yAxis = d3.svg.axis()
          .scale(x)
          .orient("left");
        
        var xAxis2 = d3.svg.axis()
          .scale(x)
          .orient("bottom");
                                
        var yAxis2 = d3.svg.axis()
          .scale(x)
          .orient("right");

        svg.selectAll ("g").data (matrix).
            enter().append ("g").selectAll("rect").data(function(d){return d;}).enter().append("rect").attr("class", "exch_rate")
            .attr("width", cellSize)
            .attr("height", cellSize)
            .attr("x", function(d,i,j) { return i * cellSize; })
            .attr("y", function(d,i,j) { return j * cellSize; })
            .style("fill", function (d) {return numberToColor(d);})
            .append("title").text(function(d,i,j) { if (typeof(d) == 'number') return d.toPrecision(7) +" (" + amino_acid_ordering[i] + "-" + amino_acid_ordering[j]+ ")"; return d;})  
                    
         svg.append("g")
          .attr("class", "x axis")
          .attr("transform", "translate(0,0)")
          .call(xAxis.tickSize(0, 0, 0));
          
      svg.append("g")
          .attr("class", "y axis")
          .attr("transform", "translate(0,0)")
          .call(yAxis.tickSize(0, 0, 0));

      svg.append("g")
          .attr("class", "x axis")
          .attr("transform", "translate(0," + 20*cellSize +")")
          .call(xAxis2.tickSize(0, 0, 0));
                
      svg.append("g")
          .attr("class", "y axis")
          .attr("transform", "translate(" + 20*cellSize +",0)")
          .call(yAxis2.tickSize(0, 0, 0));

        var where = 0;
        for (k = 0; k < aa_binning.length-1; k++) {
            where += aa_binning[k].length;
            svg.append ("line").attr ("x1",where*cellSize).attr("y1",0).attr("x2",where*cellSize).attr("y2",20*cellSize).
                                attr ('class','sep_line');
            svg.append ("line").attr ("x1",0).attr("y1",where*cellSize).attr("x2",20*cellSize).attr("y2",where*cellSize).
                                attr ('class','sep_line');
        }

        return 0;
}     

function numberToColor (value) {
    if (typeof (value) == "string") {
        return "rgba (1,0,0,1)";
    }
    rate_value = Math.min(20,Math.max(-20,Math.log (value)/Math.LN20));
    if (rate_value < 0) {
        return _colorizerB ((20+rate_value)/20.);
    }
    return _colorizerR ((rate_value)/20.);
}

function fgColor (value, normalizer) {
    if (typeof (value) == "string") {
        return "rgba (1,1,1,1)";
    }
    
    var score = 1-Math.exp(-20*value/normalizer);
    
    if (score > 0.45) {
        return "white";
    }
    return "black";
}