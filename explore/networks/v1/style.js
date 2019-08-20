vizmap = [

   {selector: "node", css: {
      "shape": "ellipse",
      "text-valign":"center",
      "text-halign":"center",
      "content": "data(id)",
      "border-width": "3px",
      "background-color": "#FFFFFF",
      "border-color":"black",
      "width": "80px",
      "height": "80px",
      "font-size":"18px"}},

   {selector:"node:selected", css: {
       "text-valign":"center",
       "text-halign":"center",
       "border-color": "black",
       "content": "data(id)",
       "border-width": "3px",
       "overlay-opacity": 0.2,
       "overlay-color": "gray"
       }},


   {selector:"edge", css: {
      "line-color": "black",
      'target-arrow-color': "black",
      'target-arrow-shape': 'triangle',
      "width": "1px",
      'curve-style': 'bezier',
      'haystack-radius': 0.5
       }},


    {"selector": "edge:selected", css: {
       "overlay-opacity": 0.2,
       "overlay-color": "gray",
       "width": "2px",
       }},


    {"selector":"edge[edgeType='repressor']", css: {
       "line-style":"solid",
       "line-color": "red",
       "width":"1px",
       "curve-style": "bezier",
       "target-arrow-shape": "tee",
       "target-arrow-color": "red"
       }},

/************
    {"selector":"edge[edgeType='activator']",
       "line-style":"solid",
       "line-color": "green",
       "width":"3px",
       "curve-style": "bezier",
       "target-arrow-shape": "arrow"
       }},

********/

   ];
