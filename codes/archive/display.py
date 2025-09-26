# matplotlib: plot the voxels near the triangular mesh
voxel_2 = voxel[voxel_for_each_vertex,:]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(voxel_2[:, 0], voxel_2[:, 1], voxel_2[:, 2], color='blue', s=1, marker='.')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_box_aspect([1,1,1])
codes.set_axes_equal.execute(ax)
plt.show()

# 