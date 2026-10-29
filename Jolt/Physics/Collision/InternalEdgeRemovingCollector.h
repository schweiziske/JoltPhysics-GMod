// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2024 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Core/QuickSort.h>
#include <Jolt/Core/STLLocalAllocator.h>
#include <Jolt/Physics/Collision/CollisionDispatch.h>

//#define JPH_INTERNAL_EDGE_REMOVING_COLLECTOR_DEBUG

#ifdef JPH_INTERNAL_EDGE_REMOVING_COLLECTOR_DEBUG
#include <Jolt/Renderer/DebugRenderer.h>
#endif // JPH_INTERNAL_EDGE_REMOVING_COLLECTOR_DEBUG

JPH_NAMESPACE_BEGIN

/// Removes internal edges from collision results. Can be used to filter out 'ghost collisions'.
/// Based on: Contact generation for meshes - Pierre Terdiman (https://www.codercorner.com/MeshContacts.pdf)
///
/// Note that this class requires that CollideSettingsBase::mActiveEdgeMode == EActiveEdgeMode::CollideWithAll
/// and CollideSettingsBase::mCollectFacesMode == ECollectFacesMode::CollectFaces.
class InternalEdgeRemovingCollector : public CollideShapeCollector
{
	static constexpr uint cMaxLocalDelayedResults = 32;
	static constexpr uint cMaxLocalVoidedFeatures = 128;
	static constexpr uint cMaxLocalSupportingFaces = 8;

	static bool				sIsPointInsideFace(Vec3Arg inPoint, const CollideShapeResult::Face &inFace, Vec3Arg inNormal, float inTolerance)
	{
		if (Square(inNormal.Dot(inPoint - inFace[0])) > Square(inTolerance))
			return false;
		for (uint j = 0; j < inFace.size(); ++j)
		{
			Vec3 edge = inFace[(j + 1) % inFace.size()] - inFace[j];
			if (edge.Cross(inPoint - inFace[j]).Dot(inNormal) < -inTolerance * edge.Length())
				return false;
		}
		return true;
	}

	/// Check if a vertex is voided
	inline bool				IsVoided(const SubShapeID &inSubShapeID, Vec3 inV) const
	{
		for (const Voided &vf : mVoidedFeatures)
			if (vf.mSubShapeID == inSubShapeID
				&& inV.IsClose(Vec3::sLoadFloat3Unsafe(vf.mFeature), mVertexToleranceSq))
				return true;
		return false;
	}

	/// Add all vertices of a face to the voided features
	inline void				VoidFeatures(const CollideShapeResult &inResult)
	{
		for (const Vec3 &v : inResult.mShape2Face)
			if (!IsVoided(inResult.mSubShapeID1, v))
			{
				Voided vf;
				v.StoreFloat3(&vf.mFeature);
				vf.mSubShapeID = inResult.mSubShapeID1;
				mVoidedFeatures.push_back(vf);
			}
	}

	/// Call the chained collector
	inline void				Chain(const CollideShapeResult &inResult)
	{
		// Make sure the chained collector has the same context as we do
		mChainedCollector.SetContext(GetContext());

		// Forward the hit
		mChainedCollector.AddHit(inResult);

		// If our chained collector updated its early out fraction, we need to follow
		UpdateEarlyOutFraction(mChainedCollector.GetEarlyOutFraction());
	}

	/// Call the chained collector and void all features of inResult
	inline void				ChainAndVoid(const CollideShapeResult &inResult)
	{
		Chain(inResult);
		VoidFeatures(inResult);
		if (mFaceTolerance > 0.0f && inResult.mShape2Face.size() >= 3)
		{
			Vec3 normal = (inResult.mShape2Face[1] - inResult.mShape2Face[0]).Cross(inResult.mShape2Face[2] - inResult.mShape2Face[0]);
			float length = normal.Length();
			if (length >= 1.0e-6f)
				mFaces.push_back({ inResult.mShape2Face, normal / length, inResult.mSubShapeID1, inResult.mSubShapeID2 });
		}

	#ifdef JPH_INTERNAL_EDGE_REMOVING_COLLECTOR_DEBUG
		DebugRenderer::sInstance->DrawWirePolygon(RMat44::sTranslation(mBaseOffset), inResult.mShape2Face, Color::sGreen);
		DebugRenderer::sInstance->DrawArrow(mBaseOffset + inResult.mContactPointOn2, mBaseOffset + inResult.mContactPointOn2 + inResult.mPenetrationAxis.NormalizedOr(Vec3::sZero()), Color::sGreen, 0.1f);
	#endif // JPH_INTERNAL_EDGE_REMOVING_COLLECTOR_DEBUG
	}

public:
	/// True for an ordinary face contact, which must not be removed as a seam.
	static bool				sIsFaceContact(const CollideShapeResult &inResult)
	{
		if (inResult.mShape2Face.size() < 3)
			return false;
		Vec3 normal = (inResult.mShape2Face[1] - inResult.mShape2Face[0]).Cross(inResult.mShape2Face[2] - inResult.mShape2Face[0]);
		float length = normal.Length();
		return length >= 1.0e-6f && normal.Dot(-inResult.mPenetrationAxis) > 0.999848f * length * inResult.mPenetrationAxis.Length();
	}

	/// Test an edge contact against a supporting face at the same moving shape pose.
	/// Also used by CCD, whose tangential casts need an overlap query to find support.
	/// The support's separating axis may itself be an edge normal: use its polygon,
	/// not that axis, to determine whether the candidate point is covered.
	static bool				sIsContactCoveredByFace(const CollideShapeResult &inEdge, const CollideShapeResult &inSupport, float inTolerance)
	{
		if (inTolerance <= 0.0f || sIsFaceContact(inEdge) || inEdge.mSubShapeID1 != inSupport.mSubShapeID1
			|| inEdge.mSubShapeID2 == inSupport.mSubShapeID2 || inSupport.mShape2Face.size() < 3)
			return false;
		const auto &face = inSupport.mShape2Face;
		Vec3 normal = (face[1] - face[0]).Cross(face[2] - face[0]);
		float length = normal.Length();
		if (length < 1.0e-6f)
			return false;
		normal /= length;
		return sIsPointInsideFace(inEdge.mContactPointOn2, face, normal, inTolerance);
	}

	/// Constructor, configures a collector to be called with all the results that do not hit internal edges
	explicit				InternalEdgeRemovingCollector(CollideShapeCollector &inChainedCollector,
			float inVertexToleranceSq
		#ifdef JPH_INTERNAL_EDGE_REMOVING_COLLECTOR_DEBUG
			, RVec3Arg inBaseOffset
		#endif // JPH_INTERNAL_EDGE_REMOVING_COLLECTOR_DEBUG
			, float inFaceTolerance = 0.0f
		) :
		CollideShapeCollector(inChainedCollector),
		mChainedCollector(inChainedCollector),
		mVertexToleranceSq(inVertexToleranceSq),
		mFaceTolerance(inFaceTolerance)
		#ifdef JPH_INTERNAL_EDGE_REMOVING_COLLECTOR_DEBUG
			, mBaseOffset(inBaseOffset)
		#endif // JPH_INTERNAL_EDGE_REMOVING_COLLECTOR_DEBUG
	{
		// Initialize arrays to full capacity to avoid needless reallocation calls
		mVoidedFeatures.reserve(cMaxLocalVoidedFeatures);
		mDelayedResults.reserve(cMaxLocalDelayedResults);
	}

	// See: CollideShapeCollector::Reset
	virtual void			Reset() override
	{
		CollideShapeCollector::Reset();

		mChainedCollector.Reset();

		mVoidedFeatures.clear();
		mDelayedResults.clear();
		mFaces.clear();
	}

	// See: CollideShapeCollector::OnBody
	virtual void			OnBody(const Body &inBody) override
	{
		// Just forward the call to our chained collector
		mChainedCollector.OnBody(inBody);
	}

	// See: CollideShapeCollector::AddHit
	virtual void			AddHit(const CollideShapeResult &inResult) override
	{
		// We only support welding when the shape is a triangle or has more vertices so that we can calculate a normal
		if (inResult.mShape2Face.size() < 3)
			return ChainAndVoid(inResult);

		// Get the triangle normal of shape 2 face
		Vec3 triangle_normal = (inResult.mShape2Face[1] - inResult.mShape2Face[0]).Cross(inResult.mShape2Face[2] - inResult.mShape2Face[0]);
		float triangle_normal_len = triangle_normal.Length();
		if (triangle_normal_len < 1e-6f)
			return ChainAndVoid(inResult);

		// If the triangle normal matches the contact normal within 1 degree, we can process the contact immediately
		// We make the assumption here that if the contact normal and the triangle normal align that the we're dealing with a 'face contact'
		Vec3 contact_normal = -inResult.mPenetrationAxis;
		float contact_normal_len = inResult.mPenetrationAxis.Length();
		if (triangle_normal.Dot(contact_normal) > 0.999848f * contact_normal_len * triangle_normal_len) // cos(1 degree)
			return ChainAndVoid(inResult);

		// Delayed processing
		mDelayedResults.push_back(inResult);
	}

	/// After all hits have been added, call this function to process the delayed results
	void					Flush()
	{
		// Sort on biggest penetration depth first
		Array<uint, STLLocalAllocator<uint, cMaxLocalDelayedResults>> sorted_indices;
		sorted_indices.resize(mDelayedResults.size());
		for (uint i = 0; i < uint(mDelayedResults.size()); ++i)
			sorted_indices[i] = i;
		QuickSort(sorted_indices.begin(), sorted_indices.end(), [this](uint inLHS, uint inRHS) { return mDelayedResults[inLHS].mPenetrationDepth > mDelayedResults[inRHS].mPenetrationDepth; });

		// Loop over all results
		for (uint i = 0; i < uint(mDelayedResults.size()); ++i)
		{
			const CollideShapeResult &r = mDelayedResults[sorted_indices[i]];

			// Determine which vertex or which edge is the closest to the contact point
			float best_dist_sq = FLT_MAX;
			uint best_v1_idx = 0;
			uint best_v2_idx = 0;
			uint num_v = uint(r.mShape2Face.size());
			uint v1_idx = num_v - 1;
			Vec3 v1 = r.mShape2Face[v1_idx] - r.mContactPointOn2;
			for (uint v2_idx = 0; v2_idx < num_v; ++v2_idx)
			{
				Vec3 v2 = r.mShape2Face[v2_idx] - r.mContactPointOn2;
				Vec3 v1_v2 = v2 - v1;
				float denominator = v1_v2.LengthSq();
				if (denominator < Square(FLT_EPSILON))
				{
					// Degenerate, assume v1 is closest, v2 will be tested in a later iteration
					float v1_len_sq = v1.LengthSq();
					if (v1_len_sq < best_dist_sq)
					{
						best_dist_sq = v1_len_sq;
						best_v1_idx = v1_idx;
						best_v2_idx = v1_idx;
					}
				}
				else
				{
					// Taken from ClosestPoint::GetBaryCentricCoordinates
					float fraction = -v1.Dot(v1_v2) / denominator;
					if (fraction < 1.0e-6f)
					{
						// Closest lies on v1
						float v1_len_sq = v1.LengthSq();
						if (v1_len_sq < best_dist_sq)
						{
							best_dist_sq = v1_len_sq;
							best_v1_idx = v1_idx;
							best_v2_idx = v1_idx;
						}
					}
					else if (fraction < 1.0f - 1.0e-6f)
					{
						// Closest lies on the line segment v1, v2
						Vec3 closest = v1 + fraction * v1_v2;
						float closest_len_sq = closest.LengthSq();
						if (closest_len_sq < best_dist_sq)
						{
							best_dist_sq = closest_len_sq;
							best_v1_idx = v1_idx;
							best_v2_idx = v2_idx;
						}
					}
					// else closest is v2, but v2 will be tested in a later iteration
				}

				v1_idx = v2_idx;
				v1 = v2;
			}

			// Check if this vertex/edge is voided
			bool voided = IsVoided(r.mSubShapeID1, r.mShape2Face[best_v1_idx])
				&& (best_v1_idx == best_v2_idx || IsVoided(r.mSubShapeID1, r.mShape2Face[best_v2_idx]));
			// Overlapping convex hulls need not share edge vertices. A rounded edge can
			// also put the closest point just outside the neighbouring supporting face.
			// Only accepted face contacts on this moving leaf can cover such an edge;
			// independent leaves and ordinary face contacts must remain untouched.
			if (!voided)
				for (const SupportingFace &face : mFaces)
				{
					if (face.mSubShapeID1 != r.mSubShapeID1 || face.mSubShapeID2 == r.mSubShapeID2)
						continue;
					if (sIsPointInsideFace(r.mContactPointOn2, face.mVertices, face.mNormal, mFaceTolerance))
					{
						voided = true;
						break;
					}
				}

		#ifdef JPH_INTERNAL_EDGE_REMOVING_COLLECTOR_DEBUG
			Color color = voided? Color::sRed : Color::sYellow;
			DebugRenderer::sInstance->DrawText3D(mBaseOffset + r.mContactPointOn2, StringFormat("%d: %g", i, r.mPenetrationDepth), color, 0.1f);
			DebugRenderer::sInstance->DrawWirePolygon(RMat44::sTranslation(mBaseOffset), r.mShape2Face, color);
			DebugRenderer::sInstance->DrawArrow(mBaseOffset + r.mContactPointOn2, mBaseOffset + r.mContactPointOn2 + r.mPenetrationAxis.NormalizedOr(Vec3::sZero()), color, 0.1f);
			DebugRenderer::sInstance->DrawMarker(mBaseOffset + r.mShape2Face[best_v1_idx], IsVoided(r.mSubShapeID1, r.mShape2Face[best_v1_idx])? Color::sRed : Color::sYellow, 0.1f);
			DebugRenderer::sInstance->DrawMarker(mBaseOffset + r.mShape2Face[best_v2_idx], IsVoided(r.mSubShapeID1, r.mShape2Face[best_v2_idx])? Color::sRed : Color::sYellow, 0.1f);
		#endif // JPH_INTERNAL_EDGE_REMOVING_COLLECTOR_DEBUG

			// No voided features, accept the contact
			if (!voided)
				Chain(r);

			// Void the features of this face
			VoidFeatures(r);
		}

		// All delayed results have been processed
		mVoidedFeatures.clear();
		mDelayedResults.clear();
		mFaces.clear();
	}

	// See: CollideShapeCollector::OnBodyEnd
	virtual void			OnBodyEnd() override
	{
		Flush();
		mChainedCollector.OnBodyEnd();
	}

	/// Version of CollisionDispatch::sCollideShapeVsShape that removes internal edges
	static void				sCollideShapeVsShape(const Shape *inShape1, const Shape *inShape2, Vec3Arg inScale1, Vec3Arg inScale2, Mat44Arg inCenterOfMassTransform1, Mat44Arg inCenterOfMassTransform2, const SubShapeIDCreator &inSubShapeIDCreator1, const SubShapeIDCreator &inSubShapeIDCreator2, const CollideShapeSettings &inCollideShapeSettings, CollideShapeCollector &ioCollector, const ShapeFilter &inShapeFilter = { }
#ifdef JPH_INTERNAL_EDGE_REMOVING_COLLECTOR_DEBUG
		, RVec3Arg inBaseOffset = RVec3::sZero()
#endif // JPH_INTERNAL_EDGE_REMOVING_COLLECTOR_DEBUG
		)
	{
		JPH_ASSERT(inCollideShapeSettings.mActiveEdgeMode == EActiveEdgeMode::CollideWithAll); // Won't work without colliding with all edges
		JPH_ASSERT(inCollideShapeSettings.mCollectFacesMode == ECollectFacesMode::CollectFaces); // Won't work without collecting faces

		InternalEdgeRemovingCollector wrapper(ioCollector, inCollideShapeSettings.mInternalEdgeRemovalVertexToleranceSq
		#ifdef JPH_INTERNAL_EDGE_REMOVING_COLLECTOR_DEBUG
			, inBaseOffset
		#endif // JPH_INTERNAL_EDGE_REMOVING_COLLECTOR_DEBUG
			, inCollideShapeSettings.mInternalEdgeRemovalFaceTolerance
		);
		CollisionDispatch::sCollideShapeVsShape(inShape1, inShape2, inScale1, inScale2, inCenterOfMassTransform1, inCenterOfMassTransform2, inSubShapeIDCreator1, inSubShapeIDCreator2, inCollideShapeSettings, wrapper, inShapeFilter);
		wrapper.Flush();
	}

private:
	// This algorithm tests a convex shape (shape 1) against a set of polygons (shape 2).
	// This assumption doesn't hold if the shape we're testing is a compound shape, so we must also
	// store the sub shape ID and ignore voided features that belong to another sub shape ID.
	struct Voided
	{
		Float3				mFeature;				// Feature that is voided (of shape 2). Read with Vec3::sLoadFloat3Unsafe so must not be the last member.
		SubShapeID			mSubShapeID;			// Sub shape ID of the shape that is colliding against the feature (of shape 1).
	};

	struct SupportingFace
	{
		CollideShapeResult::Face mVertices;
		Vec3				mNormal;
		SubShapeID			mSubShapeID1;
		SubShapeID			mSubShapeID2;
	};

	CollideShapeCollector &	mChainedCollector;
	Array<Voided, STLLocalAllocator<Voided, cMaxLocalVoidedFeatures>> mVoidedFeatures;
	Array<CollideShapeResult, STLLocalAllocator<CollideShapeResult, cMaxLocalDelayedResults>> mDelayedResults;
	Array<SupportingFace, STLLocalAllocator<SupportingFace, cMaxLocalSupportingFaces>> mFaces;
	float					mVertexToleranceSq;		// Max squared distance to consider a vertex to be the same as another vertex, used to determine if a feature is voided
	float					mFaceTolerance;			// Rounding tolerance for covered contacts, zero disables face coverage
#ifdef JPH_INTERNAL_EDGE_REMOVING_COLLECTOR_DEBUG
	RVec3					mBaseOffset; 			// Base offset for the query, used to draw the results in the right place
#endif // JPH_INTERNAL_EDGE_REMOVING_COLLECTOR_DEBUG
};

JPH_NAMESPACE_END
