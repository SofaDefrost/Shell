void main()
{
    gl_Position = ftransform();

    vec3 colour = gl_MultiTexCoord0.xyz;

    float wireframeOn = gl_Color[3];

    // If wireframe is rendered, just used white
    if (wireframeOn == 1.0)
    {
        gl_FrontColor = gl_Color;
    }
    // otherwise render as a colour map
    else
    {
        gl_FrontColor.xyz = colour;
    }
}