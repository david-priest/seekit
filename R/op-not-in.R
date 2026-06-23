# Negated match operator, migrated from CyTOF nBass_helpers.R into seekit.
# `x %!in% y` is shorthand for `!(x %in% y)`.
`%!in%` <- Negate(`%in%`)
