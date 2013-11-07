$(function() {
    $.get("outIAPWS84",
        function(d) {
            test("test", function() {
                var data = d.split("\n");
                for (var i = 0; i < data.length; i++) {
                    if (data[i] === "") {
                        break;
                    }
                    var coln = data[i].split(" ");
                    var T = parseFloat(coln[0]);
                    var sV = specific_enthalpy_of_saturated_vapor(T) * 1e-3;
                    var sL = specific_enthalpy_of_saturated_liquid(T) * 1e-3;
                    
                    equal(coln[1], sV.toFixed(5).toString(), "Compare sV");
                    equal(coln[2], sL.toFixed(8).toString(), "Compare sL");
                }
            });
     });
});


