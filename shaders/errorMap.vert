void main()
{
    gl_Position = ftransform();

    float alpha = gl_MultiTexCoord0.s;

    float wireframeOn = gl_Color[3];

    // If wireframe is rendered, just used white
    if (wireframeOn == 1.0)
    {
        gl_FrontColor = gl_Color;
    }
    // otherwise render as a colour map
    else
    {
        float scale = 3.0;
        alpha = scale * alpha;

        if (alpha<0.0)
        {
            gl_FrontColor = vec4(-alpha, 0.0, 1.0+alpha, 1.0);
        }
        else
        {
            gl_FrontColor = vec4(0.0, alpha, 1.0-alpha, 1.0);
        }
    }
}