--Keep the talent-threshold table compact when rendering Word or HTML.
function Table(tbl)
  local heading = pandoc.utils.stringify(tbl.head)
  if heading:find('Curvature change', 1, true) and
     heading:find('Endpoint tie', 1, true) then
    tbl.colspecs = {
      {pandoc.AlignDefault, 0.20},
      {pandoc.AlignRight, 0.24},
      {pandoc.AlignRight, 0.18}
    }
  end
  return tbl
end
