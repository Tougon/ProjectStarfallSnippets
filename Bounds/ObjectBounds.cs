using System.Collections;
using System.Collections.Generic;
using UnityEngine;


/// <summary>
/// An array of points that make up a shape.
/// </summary>
[System.Serializable]
public class Shape
{
    [Tooltip("Vertex points for the walkable area.")]
    // Points that make up a shape.
    public Vector2[] points = { new Vector2(10, 10), new Vector2(-10, 10), new Vector2(-10, -10), new Vector2(10, -10) };
}

/// <summary>
/// An area used for boundaries. Allows for more varied shapes than Unity's colliders.
/// Currently used for enemy walk area bounds, as this system is more useful for our needs than NavMesh.
/// </summary>
public class ObjectBounds : MonoBehaviour
{
    // Shapes that define an area.
    [SerializeField][Tooltip("Shapes for the walkable area.")]
    public List<Shape> shapes;

    // The shapes made up by the points. These shapes are created at runtime.
    public List<Poly> polys;

    // Used to prevent the mesh from moving as the attached objects move.
    Vector3 initialPos;

    // Mesh used to display walkable area in editor.
    Mesh range;

    [SerializeField]
    // Whether or not to display shapes in the editor.
    bool showShapes = true;

    /// <summary>
    /// Draws a mesh that displays the walkable area.
    /// </summary>
    private void OnDrawGizmos()
    {
        if (!showShapes || shapes == null)
            return;

        // For each shape
        foreach(Shape s in shapes)
        {
            Vector2[] points = s.points;

            // Triangulate points for the mesh
            Triangulator tr = new Triangulator(points);
            int[] triangles = tr.Triangulate();

            // Define vertices of the mesh
            Vector3[] vertices = new Vector3[points.Length];

            for (int i = 0; i < vertices.Length; i++)
                vertices[i] = new Vector3(points[i].x, 0.0f, points[i].y);


            Mesh m = new Mesh();
            m.vertices = vertices;
            m.triangles = triangles;

            m.RecalculateNormals();
            m.RecalculateBounds();

            // Draw the mesh at a set position if the game is running
            if (Application.isPlaying)
                Gizmos.DrawWireMesh(m, initialPos);
            else
                Gizmos.DrawWireMesh(m, transform.position);
        }
    }


    public void Awake()
    {
        initialPos = transform.position;
        polys = new List<Poly>();

        // Generate poly shapes from the shape vertices.
        foreach(Shape s in shapes)
            polys.Add(new Poly(s.points, new Vector2(initialPos.x, initialPos.z)));
    }


    /// <summary>
    /// Checks all shapes in <see cref="shapes"/> to see if the given position is in bounds.
    /// </summary>
    /// <param name="pos">Position to check.</param>
    /// <returns>Whether or not position is in bounds.</returns>
    public bool IsInBounds(Vector3 pos)
    {
        foreach(Poly p in polys)
        {
            if (p.ContainsPoint(new Vector2(pos.x, pos.z)))
                return true;
        }

        return false;
    }


    /// <summary>
    /// Checks all shapes in <see cref="shapes"/> to see if the given position is in bounds. Useful for checking if anything is close to the boundary limits.
    /// </summary>
    /// <param name="pos">Position to check.</param>
    /// <param name="scale">Size of the area to check. Capped from 0 to 1.</param>
    /// <returns></returns>
    public bool IsInBounds(Vector3 pos, float scale)
    {
        foreach (Poly p in polys)
        {
            if (p.ContainsPoint(new Vector2(pos.x, pos.z), new Vector2(initialPos.x, initialPos.z), scale))
                return true;
        }

        return false;
    }


    /// <summary>
    /// Checks all shapes in <see cref="shapes"/> to see if the given distance from any shape's edge. Best for checking if anything is close to the boundary limits.
    /// </summary>
    /// <param name="pos">Position to check.</param>
    /// <param name="distance">Distance from the edge to check.</param>
    /// <returns></returns>
    public bool IsDistanceFromEdge(Vector3 pos, float distance)
    {
        foreach (Poly p in polys)
        {
            Vector2 objScale = p.Size();
            Vector2 scale = new Vector2(Mathf.Clamp(objScale.x - (distance * 2), 0, objScale.x) / objScale.x,
                Mathf.Clamp(objScale.y - (distance * 2), 0, objScale.y) / objScale.y);

            if (!p.ContainsPoint(new Vector2(pos.x, pos.z), new Vector2(initialPos.x, initialPos.z), scale))
                return true;
        }

        return false;
    }


    /// <summary>
    /// Clamps the given point to the nearest point on the edge of the nearest shape in  <see cref="shapes"/> if it is out of bounds.
    /// </summary>
    /// <param name="pos">Position to clamp.</param>
    /// <returns>Clamped position.</returns>
    public Vector3 ClampPointToBounds(Vector3 pos)
    {
        if (IsInBounds(pos))
            return pos;

        Vector3 clamped = new Vector3();
        Vector2 posXZ = new Vector2(pos.x, pos.z);
        float distance = float.PositiveInfinity;

        foreach(Poly p in polys)
        {
            foreach(Vector2 v in p.pointsWorldSpace)
            {
                float d = Vector2.Distance(posXZ, v);

                if (d <= distance)
                {
                    clamped = p.ClampPointToShape(pos);
                    distance = d;
                }
            }
        }

        return clamped;
    }


    /// <summary>
    /// Clamps the edges of all <see cref="shapes"/> to fit within the provided bounds.
    /// </summary>
    /// <param name="bounds">The greater bounds.</param>
    public void ConfineToBounds(Poly bounds)
    {
        foreach (Poly p in polys)
            p.ClampShapeToShape(bounds);
    }


    public Vector3 GetRandomPointInBounds()
    {
        Vector2 point;

        if (polys == null || polys.Count == 0)
            return Vector3.zero;

        Poly p = polys[Random.Range(0, polys.Count)];
        Vector2 size = p.Size();

        point = p.centerWorldSpace + new Vector2(Random.Range(-size.x, size.x), Random.Range(-size.y, size.y));

        return ClampPointToBounds(new Vector3(point.x, 0.0f, point.y));
    }
}


/// <summary>
/// Class for various operations on an arbitrary 2D shape.
/// </summary>
public class Poly
{
    /// <summary>
    /// The positions of the shape's points.
    /// </summary>
    public Vector2[] points { get; private set; }
    /// <summary>
    /// The positions of the shape's points in world space.
    /// </summary>
    public Vector2[] pointsWorldSpace { get; private set; }
    /// <summary>
    /// The center of the shape.
    /// </summary>
    public Vector2 center { get; private set; }
    /// <summary>
    /// The center of the shape in world space.
    /// </summary>
    public Vector2 centerWorldSpace { get; private set; }
    /// <summary>
    /// The borders of the shape.
    /// </summary>
    public Vector2[] borders { get; private set; }

    /// <summary>
    /// Constructs a poly shape.
    /// </summary>
    /// <param name="p">Points of the shape in object space</param>
    /// <param name="pos">Position of object in world space</param>
    public Poly(Vector2[] p, Vector2 pos)
    {
        points = p;
        pointsWorldSpace = new Vector2[points.Length];

        // Convert the shape to world space. Otherwise, the shape will always be at the origin, and that's unhelpful.
        for (int i = 0; i < points.Length; i++)
            pointsWorldSpace[i] = new Vector2(points[i].x + pos.x, points[i].y + pos.y);

        FindCenter();
        FindCenterWorldSpace();
        SetBorders();
    }

    /// <summary>
    /// Constructs a poly shape.
    /// </summary>
    /// <param name="p">Points of the shape in object space</param>
    /// <param name="pws">Points of the shape in world space</param>
    public Poly(Vector2[] p, Vector2[] pws)
    {
        points = p;
        pointsWorldSpace = pws;

        FindCenter();
        SetBorders();
    }


    /// <summary>
    /// Finds the <see cref="center"/> of the poly shape.
    /// </summary>
    private void FindCenter()
    {
        float xMin = float.PositiveInfinity, xMax = float.NegativeInfinity, yMin = float.PositiveInfinity, yMax = float.NegativeInfinity;

        for (int i = 0; i < points.Length; i++)
        {
            xMin = points[i].x < xMin ? points[i].x : xMin;
            xMax = points[i].x > xMax ? points[i].x : xMax;
            yMin = points[i].y < yMin ? points[i].y : yMin;
            yMax = points[i].y > yMax ? points[i].y : yMax;
        }

        center = new Vector2(xMin + xMax, yMin + yMax) / 2.0f;
    }


    /// <summary>
    /// Finds the <see cref="center"/> of the poly shape.
    /// </summary>
    private void FindCenterWorldSpace()
    {
        float xMin = float.PositiveInfinity, xMax = float.NegativeInfinity, yMin = float.PositiveInfinity, yMax = float.NegativeInfinity;

        for (int i = 0; i < pointsWorldSpace.Length; i++)
        {
            xMin = pointsWorldSpace[i].x < xMin ? pointsWorldSpace[i].x : xMin;
            xMax = pointsWorldSpace[i].x > xMax ? pointsWorldSpace[i].x : xMax;
            yMin = pointsWorldSpace[i].y < yMin ? pointsWorldSpace[i].y : yMin;
            yMax = pointsWorldSpace[i].y > yMax ? pointsWorldSpace[i].y : yMax;
        }

        centerWorldSpace = new Vector2(xMin + xMax, yMin + yMax) / 2.0f;
    }

    /// <summary>
    /// Returns the max X/Y dimensions of an arbitrary 2D shape.
    /// </summary>
    /// <returns>X/Y scale</returns>
    public Vector2 Size()
    {
        float xMin = float.PositiveInfinity, xMax = float.NegativeInfinity, yMin = float.PositiveInfinity, yMax = float.NegativeInfinity;

        for(int i=0; i < points.Length; i++)
        {
            xMin = points[i].x < xMin ? points[i].x : xMin;
            xMax = points[i].x > xMax ? points[i].x : xMax;
            yMin = points[i].y < yMin ? points[i].y : yMin;
            yMax = points[i].y > yMax ? points[i].y : yMax;
        }

        Vector2 size = new Vector2(Mathf.Abs(xMin) + Mathf.Abs(xMax), Mathf.Abs(yMin) + Mathf.Abs(yMax));

        return size;
    }


    /// <summary>
    /// Finds the extents of the shape and returns them in the order xMin, yMin, xMax, yMax.
    /// </summary>
    /// <returns>An array containing the shape extents.</returns>
    public float[] Extents()
    {
        float[] extents = { float.PositiveInfinity, float.PositiveInfinity, float.NegativeInfinity, float.NegativeInfinity };

        for (int i = 0; i < pointsWorldSpace.Length; i++)
        {
            extents[0] = pointsWorldSpace[i].x < extents[0] ? pointsWorldSpace[i].x : extents[0];
            extents[1] = pointsWorldSpace[i].y < extents[1] ? pointsWorldSpace[i].y : extents[1];
            extents[2] = pointsWorldSpace[i].x > extents[2] ? pointsWorldSpace[i].x : extents[2];
            extents[3] = pointsWorldSpace[i].y > extents[3] ? pointsWorldSpace[i].y : extents[3];
        }

        return extents;
    }


    /// <summary>
    /// Defines the <see cref="borders"/> of the shape.
    /// </summary>
    private void SetBorders()
    {
        Vector2[] result = new Vector2[8];
        float[] extents = Extents();

        result[0] = new Vector2(extents[0], extents[1]);
        result[1] = new Vector2(extents[2], extents[1]);
        result[2] = new Vector2(extents[0], extents[3]);
        result[3] = new Vector2(extents[2], extents[3]);
        result[4] = new Vector2(extents[0], 0);
        result[5] = new Vector2(extents[2], 0);
        result[6] = new Vector2(0, extents[1]);
        result[7] = new Vector2(0, extents[3]);

        borders = result;
    }


    /// <summary>
    /// Checks if the given point is contained within this shape.
    /// </summary>
    /// <param name="p">The point to check.</param>
    /// <returns></returns>
    public bool ContainsPoint(Vector2 p)
    {
        return ContainsPoint(pointsWorldSpace, p);
    }


    /// <summary>
    /// Checks if the given point is contained within this shape.
    /// </summary>
    /// <param name="polyPoints">A list of points that make up an arbitary shape.</param>
    /// <param name="p">The point to check.</param>
    /// <returns></returns>
    private bool ContainsPoint(Vector2[] polyPoints, Vector2 p)
    {
        int j = polyPoints.Length - 1; 
        bool inside = false;
        for (int i = 0; i < polyPoints.Length; j = i++)
        {
            Vector2 pi = polyPoints[i];
            Vector2 pj = polyPoints[j];

            // If this check fires an odd number of times, the position is inside. Otherwise, it is not.
            if (((pi.y <= p.y && p.y < pj.y) || (pj.y <= p.y && p.y < pi.y)) &&
                (p.x < (pj.x - pi.x) * (p.y - pi.y) / (pj.y - pi.y) + pi.x))
                inside = !inside;
        }

        return inside;
    }


    /// <summary>
    /// Checks if the given point is contained within this shape scaled down.
    /// </summary>
    /// <param name="polyPoints">A list of points that make up an arbitary shape.</param>
    /// <param name="p">The point to check.</param>
    /// <param name="offset">Transform offset in X/Z.</param>
    /// <param name="scale">Scale of the shape.</param>
    /// <returns></returns>
    public bool ContainsPoint(Vector2 p, Vector2 offset, float scale)
    {
        scale = Mathf.Clamp(scale, 0, 1);

        Vector2[] newPolyPoints = new Vector2[points.Length];
        Vector2 size = (Size() / 2.0f) * (1 - scale);

		// We need to scale the points down and then apply the offset. If we scale the world space points down, nothing will be aligned right.
        for (int i = 0; i < points.Length; i++)
            newPolyPoints[i] = (points[i] * scale) + offset;//+ size + offset;

        return ContainsPoint(newPolyPoints, p);
    }


    /// <summary>
    /// Checks if the given point is contained within this shape.
    /// </summary>
    /// <param name="polyPoints">A list of points that make up an arbitary shape.</param>
    /// <param name="p">The point to check.</param>
    /// <param name="offset">Transform offset in X/Z.</param>
    /// <param name="scale">Scale of the shape.</param>
    /// <returns></returns>
    public bool ContainsPoint(Vector2 p, Vector2 offset, Vector2 scale)
    {
        Vector2[] newPolyPoints = new Vector2[points.Length];
        //Vector2 size = (Size() / 2.0f) * (scale);

        // We need to scale the points down and then apply the offset. If we scale the world space points down, nothing will be aligned right.
        for (int i = 0; i < points.Length; i++)
            newPolyPoints[i] = (points[i] * scale) + offset;//+ size + offset;

        return ContainsPoint(newPolyPoints, p);
    }


    /// <summary>
    /// Clamps a point to the nearest edge of the shape if it is out of bounds.
    /// </summary>
    /// <param name="pos">Point to clamp.</param>
    /// <returns>Clamped position.</returns>
    public Vector3 ClampPointToShape(Vector3 pos)
    {
        if (ContainsPoint(pointsWorldSpace, new Vector2(pos.x, pos.z)))
            return pos;

        int j = pointsWorldSpace.Length - 1;
        float dist = float.PositiveInfinity;
        Vector3 result = pos;

        for (int i = 0; i < pointsWorldSpace.Length; j = i++)
        {
			// Get the start and end points of a line
            Vector2 pi = pointsWorldSpace[i];
            Vector2 pj = pointsWorldSpace[j];

            // Calculate the nearest point on this line
            Vector3 linePoint = ClosestPointToLine(pos, new Vector3(pi.x, pos.y, pi.y), new Vector3(pj.x, pos.y, pj.y));
			
			// If the distance between pos and the line point is less than our current distance, update our result
            if ((pos - linePoint).sqrMagnitude <= dist * dist)//Vector3.Distance(pos, linePoint) <= dist
            {
                dist = Vector3.Distance(pos, linePoint);
                result = ClosestPointToLine(pos, new Vector3(pi.x, pos.y, pi.y), new Vector3(pj.x, pos.y, pj.y));
            }
        }

        return result;
    }


    /// <summary>
    /// Calculates the nearest point on a line starting from start to end.
    /// </summary>
    /// <param name="point">Point to clamp to the line</param>
    /// <param name="startPoint">Start point of the line</param>
    /// <param name="endPoint">End point of the line</param>
    /// <returns></returns>
    private Vector3 ClosestPointToLine(Vector3 point, Vector3 startPoint, Vector3 endPoint)
    {
        Vector3 pointTowardStart = point - startPoint;
        Vector3 startTowardEnd = (endPoint - startPoint).normalized;

        float lengthOfLine = Vector3.Distance(startPoint, endPoint);
        float dotProduct = Vector3.Dot(startTowardEnd, pointTowardStart);

		// Closest point is the beginning
        if (dotProduct <= 0)
            return startPoint;

		// Closest point is the ending
        if (dotProduct >= lengthOfLine)
            return endPoint;

		// Calculate percent distance on the line
        Vector3 thirdVector = startTowardEnd * dotProduct;

        Vector3 closestPointOnLine = startPoint + thirdVector;

        return closestPointOnLine;
    }


    /// <summary>
    /// Clamps the entire bounds of this shape to the bounds of the given shape.
    /// </summary>
    /// <param name="other">The shape to clamp to.</param>
    public void ClampShapeToShape(Poly other)
    {
        Vector2[] extents = other.borders;
        bool sizeChanged = false;

        for(int i=0; i<points.Length; i++)
        {
			// If any end point is not contained in the target shape, clamp it to that shape.
            if (!other.ContainsPoint(pointsWorldSpace[i]))
            {
                Vector3 clamp = other.ClampPointToShape(new Vector3(pointsWorldSpace[i].x, 0, pointsWorldSpace[i].y));
                Vector2 newPos = new Vector2(clamp.x, clamp.z);
                
				// We need to update the local space points as well!
                points[i] -= (pointsWorldSpace[i] - newPos);
                pointsWorldSpace[i] = newPos;

                sizeChanged = true;
            }
        }

        // Center and borders may have changed.
        if (sizeChanged)
        {
            FindCenter();
            SetBorders();
        }
    }
}


/// <summary>
/// Class that triangulates a list of points into a 2D mesh.
/// Code is from the Unify Community Wiki
/// </summary>
public class Triangulator
{
    private List<Vector2> m_points = new List<Vector2>();

    public Triangulator(Vector2[] points)
    {
        m_points = new List<Vector2>(points);
    }

    public int[] Triangulate()
    {
        List<int> indices = new List<int>();

        int n = m_points.Count;
        if (n < 3)
            return indices.ToArray();

        int[] V = new int[n];
        if (Area() > 0)
        {
            for (int v = 0; v < n; v++)
                V[v] = v;
        }
        else
        {
            for (int v = 0; v < n; v++)
                V[v] = (n - 1) - v;
        }

        int nv = n;
        int count = 2 * nv;
        for (int v = nv - 1; nv > 2;)
        {
            if ((count--) <= 0)
                return indices.ToArray();

            int u = v;
            if (nv <= u)
                u = 0;
            v = u + 1;
            if (nv <= v)
                v = 0;
            int w = v + 1;
            if (nv <= w)
                w = 0;

            if (Snip(u, v, w, nv, V))
            {
                int a, b, c, s, t;
                a = V[u];
                b = V[v];
                c = V[w];
                indices.Add(a);
                indices.Add(b);
                indices.Add(c);
                for (s = v, t = v + 1; t < nv; s++, t++)
                    V[s] = V[t];
                nv--;
                count = 2 * nv;
            }
        }

        indices.Reverse();
        return indices.ToArray();
    }

    private float Area()
    {
        int n = m_points.Count;
        float A = 0.0f;
        for (int p = n - 1, q = 0; q < n; p = q++)
        {
            Vector2 pval = m_points[p];
            Vector2 qval = m_points[q];
            A += pval.x * qval.y - qval.x * pval.y;
        }
        return (A * 0.5f);
    }

    private bool Snip(int u, int v, int w, int n, int[] V)
    {
        int p;
        Vector2 A = m_points[V[u]];
        Vector2 B = m_points[V[v]];
        Vector2 C = m_points[V[w]];
        if (Mathf.Epsilon > (((B.x - A.x) * (C.y - A.y)) - ((B.y - A.y) * (C.x - A.x))))
            return false;
        for (p = 0; p < n; p++)
        {
            if ((p == u) || (p == v) || (p == w))
                continue;
            Vector2 P = m_points[V[p]];
            if (InsideTriangle(A, B, C, P))
                return false;
        }
        return true;
    }

    private bool InsideTriangle(Vector2 A, Vector2 B, Vector2 C, Vector2 P)
    {
        float ax, ay, bx, by, cx, cy, apx, apy, bpx, bpy, cpx, cpy;
        float cCROSSap, bCROSScp, aCROSSbp;

        ax = C.x - B.x; ay = C.y - B.y;
        bx = A.x - C.x; by = A.y - C.y;
        cx = B.x - A.x; cy = B.y - A.y;
        apx = P.x - A.x; apy = P.y - A.y;
        bpx = P.x - B.x; bpy = P.y - B.y;
        cpx = P.x - C.x; cpy = P.y - C.y;

        aCROSSbp = ax * bpy - ay * bpx;
        cCROSSap = cx * apy - cy * apx;
        bCROSScp = bx * cpy - by * cpx;

        return ((aCROSSbp >= 0.0f) && (bCROSScp >= 0.0f) && (cCROSSap >= 0.0f));
    }
}
