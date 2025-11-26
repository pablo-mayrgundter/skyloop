import {describe, expect, it} from 'vitest'
import {EARTH_RADIUS_KM, LOS_ANGELES, SAN_FRANCISCO} from '../config/simulation'
import {buildGreatCircleArc, latLonToVector3} from './Globe'

const SCALE = EARTH_RADIUS_KM * 0.001

describe('latLonToVector3', () => {
  it('maps the equator to correct radius', () => {
    const point = { lat: 0, lon: 0 }
    const vec = latLonToVector3(point)
    expect(Math.abs(vec.length() - SCALE)).toBeLessThan(1e-6)
  })

  it('maps poles to y-axis', () => {
    const north = latLonToVector3({ lat: 90, lon: 0 })
    expect(north.x).toBeCloseTo(0)
    expect(north.z).toBeCloseTo(0)
  })
})

describe('buildGreatCircleArc', () => {
  it('generates consistent number of segments', () => {
    const points = buildGreatCircleArc(LOS_ANGELES, SAN_FRANCISCO, 10)
    expect(points).toHaveLength(11)
  })
})
