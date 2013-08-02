// To do the lighting in eye coordinates, apply gl_ModelViewMatrix to
// gl_Vertex and apply gl_NormalMatrix to gl_Normal. gl_NormalMatrix
// is the inverse transpose of the upper 3x3 corner of gl_ModelViewMatrix,
// which is what's required to rotate the normal into correct (ec) position.
// If the light positions are specified (in .c code) after the view volume
// transformation, then they will be stored in eye coordinates too, and
// you can access directly as gl_LightSource[0].position, etc.

// Varying vectors will be interpolated as they're passed to the fragment
// shader.

#version 120

attribute vec3 xyzcoord;
attribute vec2 uvcoord;
attribute vec3 normal;
attribute vec3 tangent;
attribute vec3 bitangent;

varying vec3 ec_vposition;

varying vec2 ts_uvcoord;
varying vec3 ts_veye;

varying mat3 TBN;

void main()
{
  vec3 ec_vnormal = normalize(gl_NormalMatrix*normal).xyz;
  vec3 ec_vtangent = normalize(gl_NormalMatrix*tangent).xyz;
  vec3 ec_vbitangent = normalize(gl_NormalMatrix*bitangent).xyz;

  TBN = transpose(mat3(
		       ec_vnormal,
		       ec_vtangent,
		       ec_vbitangent
		       ));

  ec_vposition = (gl_ModelViewMatrix*vec4(xyzcoord,1.0)).xyz;
  ts_uvcoord = uvcoord;
  
  ts_veye = TBN * normalize(-1.0 * ec_vposition);

  // State the location of the vertex
  gl_Position = gl_ProjectionMatrix*gl_ModelViewMatrix*vec4(xyzcoord,1.0);
}
