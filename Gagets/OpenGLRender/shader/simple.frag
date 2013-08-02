varying vec3 normal,lightDir;

void main()
{
  gl_FragColor = vec4(normal.x, normal.y, normal.z, 1.0);
}
