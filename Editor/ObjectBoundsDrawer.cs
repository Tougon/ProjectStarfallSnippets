using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;

/// <summary>
/// Custom editor script for object bounds. Allows for easy redrawing of bounds.
/// </summary>
[CustomEditor(typeof(ObjectBounds))]
public class ObjectBoundsDrawer : Editor
{
    public void OnSceneGUI()
    {
        ObjectBounds guiTarget = (ObjectBounds)target;

        // If there is no shape, don't draw anything.
        if (guiTarget == null || guiTarget.shapes == null)
            return;

        // For each shape in the object bounds
        foreach(Shape s in guiTarget.shapes)
        {
            // For each point in the shape
            for(int i=0; i<s.points.Length; i++)
            {
                // Draw a GUI handle on the point
                Vector3 handleDirection = Vector3.up;
                Vector3 worldSpacePos = guiTarget.transform.position + new Vector3(s.points[i].x, 0.0f, s.points[i].y);

                Handles.color = Color.blue;
                
                EditorGUI.BeginChangeCheck();

                // New point position
                Vector3 vertexPos =
                    Handles.FreeMoveHandle(worldSpacePos, Quaternion.identity, 1.0f, Vector3.zero, Handles.DotHandleCap);

                if (EditorGUI.EndChangeCheck())
                {
                    // Allows undo commands to undo vertex placement
                    Undo.RecordObject(guiTarget, "Modify shape vertex");
                    
                    // Set the point to the handle position if the handle is moved
                    vertexPos -= guiTarget.transform.position;
                    s.points[i] = new Vector2(vertexPos.x, vertexPos.z);
                }
            }
        }
    }
}
