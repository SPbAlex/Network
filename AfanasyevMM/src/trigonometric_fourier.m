clear all;

a0 = 1;
m = 100;
s = zeros( m );
s( 1 ) = a0 / 2;
for i = 2 : 2 : m
    s( i ) = 2 / ( i - 1 ) / pi;
end
stem( s, 'o' ), grid;
waitforbuttonpress;

x = -pi : 0.01 : pi;
y = a0 / 2 + sin( 0 * x );
for i = 2 : 2 : m
    plot( x, y ), grid;
    axis( [ -pi pi -0.2 1.2 ] );
    text( 2, 0, [ 'm=' num2str( i - 1 ) ] );
    y = y + 2 * sin( ( i - 1 ) .* x ) ./ ( i - 1 ) ./ pi;
    waitforbuttonpress;
end
y2 = y - a0 / 2;
for i = 2 : 2 : m
    plot( x, y2 ), grid;
    axis( [ -pi pi -0.6 0.6 ] );
    text( 2, -0.5, [ 'm=' num2str( m ) ',d=' num2str( i - 1 ) ] );
    y2 = y2 - 2 * sin( ( i - 1 ) .* x ) ./ ( i - 1 ) ./ pi;
    waitforbuttonpress;
end
