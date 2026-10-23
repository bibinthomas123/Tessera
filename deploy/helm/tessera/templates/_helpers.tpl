{{- define "tessera.fullname" -}}
{{- printf "%s-%s" .Release.Name .Chart.Name | trunc 63 | trimSuffix "-" -}}
{{- end -}}

{{- define "tessera.labels" -}}
app.kubernetes.io/name: {{ .Chart.Name }}
app.kubernetes.io/instance: {{ .Release.Name }}
app.kubernetes.io/version: {{ .Chart.AppVersion | quote }}
app.kubernetes.io/managed-by: {{ .Release.Service }}
{{- end -}}

{{- define "tessera.selector" -}}
app.kubernetes.io/name: {{ .root.Chart.Name }}
app.kubernetes.io/instance: {{ .root.Release.Name }}
app.kubernetes.io/component: {{ .component }}
{{- end -}}

{{- define "tessera.image" -}}
{{ .Values.image.repository }}:{{ .Values.image.tag | default .Chart.AppVersion }}
{{- end -}}

{{- define "tessera.dataClaim" -}}
{{ .Values.persistence.data.existingClaim | default (printf "%s-data" (include "tessera.fullname" .)) }}
{{- end -}}

{{- define "tessera.catalogClaim" -}}
{{ .Values.persistence.catalogs.existingClaim | default (printf "%s-catalogs" (include "tessera.fullname" .)) }}
{{- end -}}

{{/* Shared container settings: env, volumes, security */}}
{{- define "tessera.podCommon" -}}
securityContext:
  {{- toYaml .Values.podSecurityContext | nindent 2 }}
volumes:
  - name: data
    persistentVolumeClaim: {claimName: {{ include "tessera.dataClaim" . }}}
  - name: catalogs
    persistentVolumeClaim: {claimName: {{ include "tessera.catalogClaim" . }}, readOnly: true}
  - name: config
    configMap: {name: {{ include "tessera.fullname" . }}}
  - name: tmp
    emptyDir: {}
{{- end -}}

{{- define "tessera.containerCommon" -}}
image: {{ include "tessera.image" . }}
imagePullPolicy: {{ .Values.image.pullPolicy }}
securityContext:
  {{- toYaml .Values.securityContext | nindent 2 }}
envFrom:
  - configMapRef: {name: {{ include "tessera.fullname" . }}-env}
  - secretRef: {name: {{ .Values.existingSecret }}}
volumeMounts:
  - {name: data, mountPath: /data}
  - {name: catalogs, mountPath: /catalogs, readOnly: true}
  - {name: config, mountPath: /etc/tessera, readOnly: true}
  - {name: tmp, mountPath: /tmp}
{{- end -}}
