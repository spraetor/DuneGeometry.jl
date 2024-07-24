var documenterSearchIndex = {"docs":
[{"location":"#DuneGeometry.jl","page":"Index","title":"DuneGeometry.jl","text":"","category":"section"},{"location":"","page":"Index","title":"Index","text":"DuneGeometry is a port of the C++ Dune module dune-geometry to the julia language.","category":"page"},{"location":"","page":"Index","title":"Index","text":"The code is developed in a github repository: spraetor/DuneGeometry.jl","category":"page"},{"location":"#Module-Index","page":"Index","title":"Module Index","text":"","category":"section"},{"location":"","page":"Index","title":"Index","text":"Modules = [DuneGeometry]\nOrder   = [:constant, :type, :function, :macro]","category":"page"},{"location":"#Detailed-API","page":"Index","title":"Detailed API","text":"","category":"section"},{"location":"","page":"Index","title":"Index","text":"Modules = [DuneGeometry]\nOrder   = [:constant, :type, :function, :macro]","category":"page"},{"location":"#DuneGeometry.AffineGeometry-Union{Tuple{C}, Tuple{S}, Tuple{T}, Tuple{ReferenceElement{T}, AbstractVector{C}}} where {T<:Real, S<:Real, C<:AbstractVector{S}}","page":"Index","title":"DuneGeometry.AffineGeometry","text":"Create affine geometry from reference element and a vector of vertex coordinates.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.AffineGeometry-Union{Tuple{C}, Tuple{T}, Tuple{GeometryType, AbstractVector{C}}} where {T<:Real, C<:AbstractVector{T}}","page":"Index","title":"DuneGeometry.AffineGeometry","text":"Create affine geometry from GeometryType and a vector of vertex coordinates.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.AffineGeometry-Union{Tuple{T}, Tuple{GeometryType, AbstractVector{T}, AbstractMatrix{T}}} where T<:Real","page":"Index","title":"DuneGeometry.AffineGeometry","text":"Create affine geometry from GeometryType, one vertex, and the Jacobian matrix.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.AffineGeometry-Union{Tuple{T}, Tuple{ReferenceElement{T}, AbstractVector{T}, AbstractMatrix{T}}} where T<:Real","page":"Index","title":"DuneGeometry.AffineGeometry","text":"Create affine geometry from reference element, one vertex, and the Jacobian matrix.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.GeometryType","page":"Index","title":"DuneGeometry.GeometryType","text":"Unique label for each type of entities that can occur in a grid.\n\n\n\n\n\n","category":"type"},{"location":"#DuneGeometry.MultiLinearGeometry-Union{Tuple{C}, Tuple{T}, Tuple{GeometryType, AbstractVector{C}}} where {T<:Real, C<:AbstractVector{T}}","page":"Index","title":"DuneGeometry.MultiLinearGeometry","text":"Create multilinear geometry from GeometryType and a vector of vertex coordinates.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.MultiLinearGeometry-Union{Tuple{C}, Tuple{T}, Tuple{ReferenceElement{T}, AbstractVector{C}}} where {T<:Real, C<:AbstractVector{T}}","page":"Index","title":"DuneGeometry.MultiLinearGeometry","text":"Create multilinear geometry from reference element and a vector of vertex coordinates.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.ReferenceElement","page":"Index","title":"DuneGeometry.ReferenceElement","text":"ReferenceElement{T}\n\nThis class provides access to geometric and topological properties of a reference element.\n\nThis includes the number of subentities, the volume, and a method for checking if a point is contained in the reference element. The embedding of each subentity into the reference element is also provided.\n\n\n\n\n\n","category":"type"},{"location":"#Base.:==-Tuple{GeometryType, GeometryType}","page":"Index","title":"Base.:==","text":"Check for equality. This method knows that in dimension 0 and 1 all BasicTypes are equal.\n\n\n\n\n\n","category":"method"},{"location":"#Base.show-Tuple{IO, GeometryType}","page":"Index","title":"Base.show","text":"Show the GeometryType on the Standard output.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.affine-Union{Tuple{AffineGeometry{T}}, Tuple{T}} where T<:Real","page":"Index","title":"DuneGeometry.affine","text":"Always true: this is an affine geometry.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.affine-Union{Tuple{MultiLinearGeometry{T}}, Tuple{T}} where T<:Real","page":"Index","title":"DuneGeometry.affine","text":"Is this mapping affine?\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.basicType-Tuple{GeometryType}","page":"Index","title":"DuneGeometry.basicType","text":"Return the BasicType associated to the Geometry\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.basicType-Tuple{String}","page":"Index","title":"DuneGeometry.basicType","text":"Convert a string into a BasicType enum value.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.center-Union{Tuple{AffineGeometry{T}}, Tuple{T}} where T<:Real","page":"Index","title":"DuneGeometry.center","text":"Obtain the centroid of the mapping's image.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.center-Union{Tuple{MultiLinearGeometry{T}}, Tuple{T}} where T<:Real","page":"Index","title":"DuneGeometry.center","text":"Obtain the centroid of the mapping's image.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.checkInside-Union{Tuple{T}, Tuple{ReferenceElement{T}, DenseVector{T}}} where T<:Real","page":"Index","title":"DuneGeometry.checkInside","text":"checkInside(r, local)\n\nCheck if a coordinate is in the reference element.\n\nThis method returns true if the given local coordinate is within this reference element.\n\nArguments\n\nx::DenseVector: Coordinates of the point.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.corner-Union{Tuple{T}, Tuple{AffineGeometry{T}, Integer}} where T<:Real","page":"Index","title":"DuneGeometry.corner","text":"Obtain coordinates of the i-th corner.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.corner-Union{Tuple{T}, Tuple{MultiLinearGeometry{T}, Integer}} where T<:Real","page":"Index","title":"DuneGeometry.corner","text":"Obtain coordinates of the i-th corner.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.corners-Union{Tuple{AffineGeometry{T}}, Tuple{T}} where T<:Real","page":"Index","title":"DuneGeometry.corners","text":"Obtain number of corners of the corresponding reference element.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.corners-Union{Tuple{MultiLinearGeometry{T}}, Tuple{T}} where T<:Real","page":"Index","title":"DuneGeometry.corners","text":"Obtain number of corners of the corresponding reference element.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.geometry-Union{Tuple{G}, Tuple{T}, Tuple{Type{G}, ReferenceElement{T}, Integer, Integer}} where {T<:Real, G<:DuneGeometry.AbstractGeometry{T}}","page":"Index","title":"DuneGeometry.geometry","text":"geometry{G}(r, i)\n\nObtain the embedding of subentity (i,codim) into the reference element\n\nDenote by E the i-th subentity of codimension codim of the current reference element. This method returns a geometry of type G that maps the reference element of E into the current reference element.\n\nArguments\n\ni::Int: number of subentity E (0 < i <= size( c ))\nc::Int: Codimension of subentity E\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.integrationOuterNormal-Union{Tuple{T}, Tuple{ReferenceElement{T}, Integer}} where T<:Real","page":"Index","title":"DuneGeometry.integrationOuterNormal","text":"integrationOuterNormal(r, face)\n\nObtain the integration outer normal of the reference element.\n\nThe integration outer normal is the outer normal whose length coincides with the face's integration element.\n\nArguments\n\nface::Int: Index of the face, whose normal is desired.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.isConical-Tuple{GeometryType, Int64}","page":"Index","title":"DuneGeometry.isConical","text":"Return true if entity was constructed with a conical product in the chosen step.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.isConical-Tuple{GeometryType}","page":"Index","title":"DuneGeometry.isConical","text":"Return true if entity was constructed with a conical product in the last step.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.isPrismatic-Tuple{GeometryType, Int64}","page":"Index","title":"DuneGeometry.isPrismatic","text":"Return true if entity was constructed with a prismatic product in the chosen step.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.isPrismatic-Tuple{GeometryType}","page":"Index","title":"DuneGeometry.isPrismatic","text":"Return true if entity was constructed with a prismatic product in the last step.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.jacobian-Union{Tuple{S}, Tuple{T}, Tuple{AffineGeometry{T}, AbstractVector{S}}} where {T<:Real, S<:Real}","page":"Index","title":"DuneGeometry.jacobian","text":"jacobian(g,x)\n\nObtain the Jacobian.\n\nArguments:\n\nx::AbstractVector{S}: local coordinate to evaluate Jacobian in\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.jacobian-Union{Tuple{S}, Tuple{T}, Tuple{MultiLinearGeometry{T}, AbstractVector{S}}} where {T<:Real, S<:Real}","page":"Index","title":"DuneGeometry.jacobian","text":"jacobian(g,x)\n\nObtain the Jacobian.\n\nArguments:\n\nx::AbstractVector{S}: local coordinate to evaluate Jacobian in\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.jacobianTransposed-Union{Tuple{S}, Tuple{T}, Tuple{AffineGeometry{T}, AbstractVector{S}}} where {T<:Real, S<:Real}","page":"Index","title":"DuneGeometry.jacobianTransposed","text":"jacobianTransposed(g,x)\n\nObtain the transposed of the Jacobian.\n\nArguments:\n\nx::AbstractVector{S}: local coordinate to evaluate Jacobian in\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.jacobianTransposed-Union{Tuple{S}, Tuple{T}, Tuple{MultiLinearGeometry{T}, AbstractVector{S}}} where {T<:Real, S<:Real}","page":"Index","title":"DuneGeometry.jacobianTransposed","text":"jacobianTransposed(g,x)\n\nObtain the transposed of the Jacobian.\n\nArguments:\n\nx::AbstractVector{S}: local coordinate to evaluate Jacobian in\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.position-Union{Tuple{T}, Tuple{ReferenceElement{T}, Integer, Integer}} where T<:Real","page":"Index","title":"DuneGeometry.position","text":"position(r, i, c)\n\nPosition of the barycenter of entity (i,c).\n\nDenote by E the i-th subentity of codimension c of the current reference element. This method returns the coordinates of the center of gravity of E within the current reference element.\n\nArguments\n\nr::ReferenceElement: The reference element.\ni::Int: Number of subentity E (0 < i <= size( c ))\nc::Int: Codimension of subentity E\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.referenceElement-Union{Tuple{AffineGeometry{T}}, Tuple{T}} where T<:Real","page":"Index","title":"DuneGeometry.referenceElement","text":"Obtain the reference element the geometry is defined on.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.referenceElement-Union{Tuple{MultiLinearGeometry{T}}, Tuple{T}} where T<:Real","page":"Index","title":"DuneGeometry.referenceElement","text":"Obtain the reference element the geometry is defined on.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.size-Union{Tuple{T}, Tuple{ReferenceElement{T}, Integer, Integer, Integer}} where T<:Real","page":"Index","title":"DuneGeometry.size","text":"size(r, i, c, cc)\n\nNumber of subentities of codimension cc of subentity (i,c).\n\nDenote by E the i-th subentity of codimension c of the given reference element r. This method returns the number of subentities of codimension cc of the reference element, that are also a subentity of E. If cc<c this number is zero.\n\nArguments\n\nr::ReferenceElement: The reference element.\ni:Int: The number of subentity E (0 < i <= size(r,c))\nc::Int: Codimension of subentity E (0 <= c <= dim(r))\ncc::Int: Codimension whose size is desired (0 <= cc <= dim(r))\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.size-Union{Tuple{T}, Tuple{ReferenceElement{T}, Integer}} where T<:Real","page":"Index","title":"DuneGeometry.size","text":"size(r, c)\n\nNumber of subentities of codimension c in the ReferenceElement r.\n\nArguments\n\nr::ReferenceElement: The reference element.\nc::Int: Codimension whose size is requested.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.subEntities-Union{Tuple{T}, Tuple{ReferenceElement{T}, Integer, Integer, Integer}} where T<:Real","page":"Index","title":"DuneGeometry.subEntities","text":"subEntities(r, i, c, cc)\n\nObtain the range of numbers of subentities with codim cc of (i,c).\n\nDenote by E the i-th subentity of codimension c of the current reference element. This method returns a range of numbers of all subentities of E with codimension cc. Notice that the sub-subentity codimension as well as the numbers in the returned range are given with respect to the reference element itself and not with respect to E. For 0<=cc<c this will return an empty range.\n\nArguments\n\nr::ReferenceElement: The reference element.\ni::Int: Number of subentity E (0 < i <= size( c ))\nc::Int: Codimension of subentity E\ncc:Int: Codimension of subentity S (c <= cc <= dim)\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.subEntity-Union{Tuple{T}, Tuple{ReferenceElement{T}, Vararg{Integer, 4}}} where T<:Real","page":"Index","title":"DuneGeometry.subEntity","text":"subEntity(r, i, c, ii, cc)\n\nObtain number of ii-th subentity with codim cc of (i,c).\n\nDenote by E the i-th subentity of codimension c of the current reference element. And denote by S the ii-th subentity of codimension (cc-c) of E. Then, S is a also a subentity of codimension cc of the current reference element. This method returns the number of S with respect to the current reference element.\n\nArguments\n\nr::ReferenceElement: The reference element.\ni::Int: Number of subentity E (0 < i <= size( c ))\nc::Int: Codimension of subentity E\nii:Int: Number of subentity S (with respect to E)\ncc:Int: Codimension of subentity S (c <= cc <= dim)\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.toGlobal-Union{Tuple{S}, Tuple{T}, Tuple{AffineGeometry{T}, AbstractVector{S}}} where {T<:Real, S<:Real}","page":"Index","title":"DuneGeometry.toGlobal","text":"toGlobal(g,x)\n\nEvaluate the local-to-global mapping.\n\nArguments:\n\nx::AbstractVector{S}: local coordinate to map\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.toGlobal-Union{Tuple{S}, Tuple{T}, Tuple{MultiLinearGeometry{T}, AbstractVector{S}}} where {T<:Real, S<:Real}","page":"Index","title":"DuneGeometry.toGlobal","text":"toGlobal(g,x)\n\nEvaluate the local-to-global mapping.\n\nArguments:\n\nx::AbstractVector{S}: local coordinate to map\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.toId-Tuple{GeometryType}","page":"Index","title":"DuneGeometry.toId","text":"Construct an Id representing this GeometryType\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.toLocal-Union{Tuple{S}, Tuple{T}, Tuple{AffineGeometry{T}, AbstractVector{S}}} where {T<:Real, S<:Real}","page":"Index","title":"DuneGeometry.toLocal","text":"toLocal(g,X)\n\nEvaluate the inverse mapping, i.e. the global-to-local mapping.\n\nThe returned local coordinate y minimizes\n\n(toGlobal( y ) - x).two_norm()\n\non the entire affine hull of the reference element.  This degenerates to the inverse map if the argument y is in the range of the map.\n\nArguments:\n\nX::AbstractVector{S}: global coordinate to map\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.toLocal-Union{Tuple{S}, Tuple{T}, Tuple{MultiLinearGeometry{T}, AbstractVector{S}}} where {T<:Real, S<:Real}","page":"Index","title":"DuneGeometry.toLocal","text":"toLocal(g,X)\n\nEvaluate the inverse mapping, i.e. the global-to-local mapping.\n\nThe returned local coordinate y minimizes\n\n(toGlobal( y ) - x).two_norm()\n\non the entire affine hull of the reference element.  This degenerates to the inverse map if the argument y is in the range of the map.\n\nArguments:\n\nX::AbstractVector{S}: global coordinate to map\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.toString-Tuple{DuneGeometry.BasicType.T}","page":"Index","title":"DuneGeometry.toString","text":"Convert a BasicType into a string.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.toString-Tuple{GeometryType}","page":"Index","title":"DuneGeometry.toString","text":"Convert a GeometryType into a string.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.type-Union{Tuple{AffineGeometry{T}}, Tuple{T}} where T<:Real","page":"Index","title":"DuneGeometry.type","text":"Obtain the type of the reference element.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.type-Union{Tuple{MultiLinearGeometry{T}}, Tuple{T}} where T<:Real","page":"Index","title":"DuneGeometry.type","text":"Obtain the type of the reference element.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.type-Union{Tuple{ReferenceElement{T}}, Tuple{T}} where T<:Real","page":"Index","title":"DuneGeometry.type","text":"Obtain the type of this reference element\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.type-Union{Tuple{T}, Tuple{ReferenceElement{T}, Integer, Integer}} where T<:Real","page":"Index","title":"DuneGeometry.type","text":"type(r, i, c)\n\nObtain the type of subentity (i,c).\n\nDenote by E the i-th subentity of codimension c of the current reference element. This method returns the GeometryType of E.\n\nArguments\n\nr::ReferenceElement: The reference element.\ni::Int: Number of subentity E (0 < i <= size( c ))\nc::Int: Codimension of subentity E\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.volume-Union{Tuple{AffineGeometry{T}}, Tuple{T}} where T<:Real","page":"Index","title":"DuneGeometry.volume","text":"Obtain the volume of the element.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.volume-Union{Tuple{MultiLinearGeometry{T}}, Tuple{T}} where T<:Real","page":"Index","title":"DuneGeometry.volume","text":"Obtain the volume of the element.\n\n\n\n\n\n","category":"method"},{"location":"#DuneGeometry.volume-Union{Tuple{ReferenceElement{T}}, Tuple{T}} where T<:Real","page":"Index","title":"DuneGeometry.volume","text":"Obtain the volume of the reference element.\n\n\n\n\n\n","category":"method"}]
}