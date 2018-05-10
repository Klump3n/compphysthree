(TeX-add-style-hook
 "cp3_ue_2"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "10pt" "a4paper")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "a4paper" "top=2cm" "bottom=2.5cm" "left=3cm" "right=3cm") ("babel" "english" "ngerman") ("inputenc" "utf8") ("fontenc" "T1") ("newtxmath" "libertine")))
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "geometry"
    "babel"
    "inputenc"
    "fontenc"
    "graphicx"
    "lmodern"
    "epstopdf"
    "caption"
    "verbatim"
    "multicol"
    "hyphsubst"
    "xcolor"
    "soul"
    "float"
    "array"
    "subfigure"
    "ulem"
    "siunitx"
    "hyperref"
    "booktabs"
    "etex"
    "pgfplots"
    "tikz"
    "tikz-3dplot"
    "tikzscale"
    "listings"
    "amsfonts"
    "amsmath"
    "amsthm"
    "amssymb"
    "cancel"
    "mathcomp"
    "nicefrac"
    "libertine"
    "newtxmath")
   (LaTeX-add-lengths
    "figureheight"
    "figurewidth")
   (LaTeX-add-xcolor-definecolors
    "colKeys"
    "colIdentifier"
    "colComments"
    "dkgreen"
    "gray"
    "colString")
   (LaTeX-add-array-newcolumntypes
    "C"))
 :latex)

