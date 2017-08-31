PRO abline, value, $
    VERTICAL=vertical, $
    LINESTYLE=linestyle, $
    THICK=thick, $
    COLOR=color
    
   IF N_Elements(color) EQ 0 THEN color = 'white'

   IF Keyword_Set(Vertical) THEN BEGIN
       PLOTS, [value, value], [!Y.CRange[0], !Y.CRange[1]], $
          LINESTYLE=linestyle, THICK=thick, $
          COLOR=FSC_Color(String(color))   
   ENDIF ELSE BEGIN
       PLOTS, [!X.CRange[0], !X.CRange[1]], [value, value], $
          LINESTYLE=linestyle, THICK=thick, $
          COLOR=FSC_Color(color)
   ENDELSE

END