package mtGeometry

import (
	"math"
	"unsafe"

	"github.com/go-gl/gl/v3.3-core/gl"
	"github.com/go-gl/mathgl/mgl32"
	//"fmt"
)

type Geometry struct {
	// Important geometry attributes
	ArrayBuffer  uint32
	VertexBuffer uint32
	IndexBuffer  uint32
	IndexCount   int32
}

type ArrayGeometry struct {
	ArrayBuffer  uint32
	VertexBuffer uint32
	VertexCount  int32
}

type Object struct {
	Geo     Geometry
	Pos     mgl32.Vec3
	Scale   mgl32.Vec3
	Color   mgl32.Vec3
	IsLight bool
}

type Mesh struct {
	Pos    mgl32.Vec3
	Normal mgl32.Vec3
	UV     mgl32.Vec2
}

func createUnitSphere(subdivLat, subdivLong int, mesh *[]Mesh, meshIndices *[]uint32) {
	var radius float32 = 0.5

	// +1.0f because there's a gap between the poles and the first parallel.
	latitudeSpacing := 1.0 / (float32(subdivLat) + 1.0)
	longitudeSpacing := 1.0 / float32(subdivLong)

	v := 0
	// North pole.
	for longitude := 0; longitude <= subdivLong+1; longitude++ {
		(*mesh)[v].Pos = mgl32.Vec3{0, radius, 0}
		(*mesh)[v].Normal = mgl32.Vec3{0, radius, 0}.Normalize()
		(*mesh)[v].UV = mgl32.Vec2{float32(longitude) * longitudeSpacing, 1}
		v++
	}

	for latitude := 0; latitude < subdivLat; latitude++ {
		for longitude := 0; longitude <= subdivLong; longitude++ {

			// Scale coordinates into the 0...1 texture coordinate range,
			// with north at the top (y = 1).
			(*mesh)[v].UV = mgl32.Vec2{float32(longitude) * longitudeSpacing, 1.0 - float32(latitude+1)*latitudeSpacing}

			// Convert to spherical coordinates:
			// theta is a longitude angle (around the equator) in radians.
			// phi is a latitude angle (north or south of the equator).
			theta := float64((*mesh)[v].UV.X()) * 2.0 * math.Pi
			phi := (float64((*mesh)[v].UV.Y()) - 0.5) * math.Pi

			// This determines the radius of the ring of this line of latitude.
			// It's widest at the equator, and narrows as phi increases/decreases.
			c := float32(math.Cos(phi))

			// Usual formula for a vector in spherical coordinates.
			// You can exchange x & z to wind the opposite way around the sphere.
			(*mesh)[v].Pos = mgl32.Vec3{c * float32(math.Cos(theta)), float32(math.Sin(phi)), c * float32(math.Sin(theta))}.Mul(radius)
			(*mesh)[v].Normal = (*mesh)[v].Pos.Normalize()

			v++
		}
	}

	// South pole.
	for longitude := 0; longitude <= subdivLong+1; longitude++ {
		(*mesh)[v].Pos = mgl32.Vec3{0, -radius, 0}
		(*mesh)[v].Normal = mgl32.Vec3{0, -radius, 0}.Normalize()
		(*mesh)[v].UV = mgl32.Vec2{float32(longitude) * longitudeSpacing, 0}
		v++
	}

	i := 0
	// Triangles for the mesh in between the poles. This should be just like any normal rectangular mesh.
	for latitude := uint32(0); latitude <= uint32(subdivLat); latitude++ {
		for longitude := uint32(0); longitude < uint32(subdivLong); longitude++ {
			(*meshIndices)[i] = 1 + (latitude * (uint32(subdivLong) + 1)) + longitude
			(*meshIndices)[i+1] = 1 + (latitude * (uint32(subdivLong) + 1)) + longitude + uint32(subdivLong) + 1
			(*meshIndices)[i+2] = 1 + (latitude * (uint32(subdivLong) + 1)) + longitude + 1

			(*meshIndices)[i+3] = 1 + (latitude * (uint32(subdivLong) + 1)) + longitude + 1
			(*meshIndices)[i+4] = 1 + (latitude * (uint32(subdivLong) + 1)) + longitude + uint32(subdivLong) + 1
			(*meshIndices)[i+5] = 1 + (latitude * (uint32(subdivLong) + 1)) + longitude + uint32(subdivLong) + 2
			i += 6
		}
	}
}

// Defines the pure vertex and index arrays. Given the pre-allocated arrays, it directly writes into them at
// a given position.
func createUnitSquare(subdivisions int, offset mgl32.Vec3, mesh *[]Mesh, meshOffset int, meshIndices *[]uint32, meshIndexOffset int) {

	size := subdivisions + 2

	sizeFactor := 1.0 / float32(size-1)

	for i := 0; i < size; i++ {
		for j := 0; j < size; j++ {
			// It's now centered at the origin on both axes and has unit size.
			// Possibly with an offset defined.
			(*mesh)[j*size+i].Pos = mgl32.Vec3{float32(i)*sizeFactor - 0.5 + offset.X(), offset.Y(), float32(j)*sizeFactor - 0.5 + offset.Z()}
			(*mesh)[j*size+i].Normal = mgl32.Vec3{0, 1, 0}
			(*mesh)[j*size+i].UV = mgl32.Vec2{float32(i) * sizeFactor, float32(j) * sizeFactor}
		}
	}

	index := 0
	uSize := uint32(size)
	for i := uint32(0); i < uSize-1; i++ {
		for j := uint32(0); j < uSize-1; j++ {
			(*meshIndices)[index] = j*uSize + i
			(*meshIndices)[index+1] = j*uSize + i + 1
			(*meshIndices)[index+2] = (j+1)*uSize + i + 1

			(*meshIndices)[index+3] = j*uSize + i
			(*meshIndices)[index+4] = (j+1)*uSize + i + 1
			(*meshIndices)[index+5] = (j+1)*uSize + i
			index += 6
		}
	}

}

func generateGeometryAttributes(mesh *[]Mesh, meshIndices *[]uint32, vertexCount, indexCount int) Geometry {
	geo := Geometry{}

	var m Mesh
	stride := int32(unsafe.Sizeof(m))
	var v mgl32.Vec3
	vStride := int(unsafe.Sizeof(v))
	var ui uint32
	uiStride := int(unsafe.Sizeof(ui))

	gl.GenBuffers(1, &geo.ArrayBuffer)
	gl.BindBuffer(gl.ARRAY_BUFFER, geo.ArrayBuffer)
	gl.BufferData(gl.ARRAY_BUFFER, int(stride)*vertexCount, gl.Ptr(*mesh), gl.STATIC_DRAW)

	gl.GenVertexArrays(1, &geo.VertexBuffer)
	gl.BindVertexArray(geo.VertexBuffer)

	gl.EnableVertexAttribArray(0)
	gl.VertexAttribPointer(0, 3, gl.FLOAT, false, stride, gl.PtrOffset(0))
	gl.EnableVertexAttribArray(1)
	gl.VertexAttribPointer(1, 3, gl.FLOAT, true, stride, gl.PtrOffset(vStride))
	gl.EnableVertexAttribArray(2)
	gl.VertexAttribPointer(2, 2, gl.FLOAT, false, stride, gl.PtrOffset(2*vStride))

	gl.BindBuffer(gl.ARRAY_BUFFER, 0)

	gl.GenBuffers(1, &geo.IndexBuffer)
	gl.BindBuffer(gl.ELEMENT_ARRAY_BUFFER, geo.IndexBuffer)
	gl.BufferData(gl.ELEMENT_ARRAY_BUFFER, uiStride*indexCount, gl.Ptr(*meshIndices), gl.STATIC_DRAW)
	gl.BindBuffer(gl.ELEMENT_ARRAY_BUFFER, 0)

	geo.IndexCount = int32(indexCount)

	return geo
}

func GenerateGeometryAttributes(mesh *[]Mesh, meshIndices *[]uint32, vertexCount, indexCount int) Geometry {
	return generateGeometryAttributes(mesh, meshIndices, vertexCount, indexCount)
}

func GenerateGeometryArrayAttributes(mesh *[]Mesh, vertexCount int) ArrayGeometry {
	geo := ArrayGeometry{}

	var m Mesh
	stride := int32(unsafe.Sizeof(m))
	var v mgl32.Vec3
	vStride := int(unsafe.Sizeof(v))

	gl.GenBuffers(1, &geo.ArrayBuffer)
	gl.BindBuffer(gl.ARRAY_BUFFER, geo.ArrayBuffer)
	gl.BufferData(gl.ARRAY_BUFFER, int(stride)*vertexCount, gl.Ptr(*mesh), gl.STATIC_DRAW)

	gl.GenVertexArrays(1, &geo.VertexBuffer)
	gl.BindVertexArray(geo.VertexBuffer)

	gl.EnableVertexAttribArray(0)
	gl.VertexAttribPointer(0, 3, gl.FLOAT, false, stride, gl.PtrOffset(0))
	gl.EnableVertexAttribArray(1)
	gl.VertexAttribPointer(1, 3, gl.FLOAT, true, stride, gl.PtrOffset(vStride))
	gl.EnableVertexAttribArray(2)
	gl.VertexAttribPointer(2, 2, gl.FLOAT, false, stride, gl.PtrOffset(2*vStride))

	gl.BindBuffer(gl.ARRAY_BUFFER, 0)

	geo.VertexCount = int32(vertexCount)

	return geo
}

func CreateFullscreenQuadGeometry() Geometry {
	mesh := []Mesh{
		Mesh{mgl32.Vec3{-1.0, -1.0, 0.0}, mgl32.Vec3{0.0, 0.0, 1.0}, mgl32.Vec2{0.0, 0.0}},
		Mesh{mgl32.Vec3{1.0, -1.0, 0.0}, mgl32.Vec3{0.0, 0.0, 1.0}, mgl32.Vec2{1.0, 0.0}},
		Mesh{mgl32.Vec3{1.0, 1.0, 0.0}, mgl32.Vec3{0.0, 0.0, 1.0}, mgl32.Vec2{1.0, 1.0}},
		Mesh{mgl32.Vec3{-1.0, 1.0, 0.0}, mgl32.Vec3{0.0, 0.0, 1.0}, mgl32.Vec2{0.0, 1.0}},
	}
	meshIndices := []uint32{0, 1, 2, 0, 2, 3}

	return generateGeometryAttributes(&mesh, &meshIndices, 4, 6)
}

func CreateUnitSphereGeometry(subdivLat, subdivLong int) Geometry {
	vertexCount := (subdivLat+2)*(subdivLong+1) + 2
	mesh := make([]Mesh, vertexCount, vertexCount)

	// (subdivLong+1)*3*2 = A triangle from the first vertex to each of the first latitute.
	// for both poles.
	indexCount := (subdivLat + 1) * subdivLong * 6
	meshIndices := make([]uint32, indexCount, indexCount)

	createUnitSphere(subdivLat, subdivLong, &mesh, &meshIndices)

	return generateGeometryAttributes(&mesh, &meshIndices, vertexCount, indexCount)

}

func CreateUnitSquareGeometry(subdivisions int, offset mgl32.Vec3) Geometry {
	size := subdivisions + 2
	vertexCount := size * size
	mesh := make([]Mesh, vertexCount, vertexCount)

	indexCount := (size - 1) * (size - 1) * 2 * 3
	meshIndices := make([]uint32, indexCount, indexCount)

	createUnitSquare(subdivisions, offset, &mesh, 0, &meshIndices, 0)

	return generateGeometryAttributes(&mesh, &meshIndices, vertexCount, indexCount)
}

func CreateObject(geo Geometry, pos, scale, color mgl32.Vec3, isLight bool) Object {
	return Object{
		Geo:     geo,
		Pos:     pos,
		Scale:   scale.Mul(0.5),
		Color:   color,
		IsLight: isLight,
	}
}
